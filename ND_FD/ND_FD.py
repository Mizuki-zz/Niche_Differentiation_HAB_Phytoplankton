import pandas as pd
import numpy as np
import os


# Define calculation function
def calculate_coexistence_metrics(beta_matrix, abundance_df):
    """
    Calculate niche differences, fitness differences, and coexistence metrics

    Parameters:
    beta_matrix: species interaction matrix (DataFrame)
    abundance_df: species abundance data (DataFrame)

    Returns:
    DataFrame containing Ni, Fi, and coexistence judgment
    """
    species = beta_matrix.columns.tolist()
    n_species = len(species)

    # Calculate maximum abundance N_j^* for each species
    max_abundance = abundance_df.max()

    # Calculate when species i is at its minimum, the abundance N_j^{-i,*} of species j
    min_abundance_matrix = pd.DataFrame(index=species, columns=species)

    for i in species:
        # Find the sample with minimum abundance of species i
        min_idx = abundance_df[i].idxmin()
        # Record the abundance of all species in that sample
        min_abundance_matrix.loc[i] = abundance_df.loc[min_idx]

    # Initialize result DataFrame
    results = pd.DataFrame(index=species, columns=['Ni', 'Fi', 'Coexistence'])

    # Calculate nij and fij for each species pair
    nij_matrix = pd.DataFrame(np.zeros((n_species, n_species)), index=species, columns=species)
    fij_matrix = pd.DataFrame(np.zeros((n_species, n_species)), index=species, columns=species)
    cij_matrix = pd.DataFrame(np.zeros((n_species, n_species)), index=species, columns=species)

    for i in species:
        for j in species:
            if i == j:
                continue

            beta_ii = beta_matrix.loc[i, i]
            beta_jj = beta_matrix.loc[j, j]
            beta_ij = beta_matrix.loc[i, j]
            beta_ji = beta_matrix.loc[j, i]

            # Calculate nij (Equation 8)
            numerator_nij = np.abs(beta_ij * beta_ji)
            denominator_nij = np.abs(beta_jj * beta_ii)
            if denominator_nij > 0:
                nij = 1 - np.sqrt(numerator_nij / denominator_nij)
            else:
                nij = np.nan
            nij_matrix.loc[i, j] = nij

            # Calculate fij (Equation 9)
            numerator_fij = np.abs(beta_ji * beta_jj)
            denominator_fij = np.abs(beta_ij * beta_ii)
            if denominator_fij > 0:
                fij = np.sqrt(numerator_fij / denominator_fij)
            else:
                fij = np.nan
            fij_matrix.loc[i, j] = fij

            # Calculate cij (Equation 12)
            numerator_cij = np.abs(beta_ii * beta_ji)
            denominator_cij = np.abs(beta_jj * beta_ij)
            if denominator_cij > 0:
                cij = np.sqrt(numerator_cij / denominator_cij)
            else:
                cij = np.nan
            cij_matrix.loc[i, j] = cij

    # Calculate Ni and Fi for each species
    for i in species:
        numerator_ni = 0
        denominator_ni = 0
        fi_sum = 0

        for j in species:
            if i == j:
                continue

            cij = cij_matrix.loc[i, j]
            nij = nij_matrix.loc[i, j]
            fij = fij_matrix.loc[i, j]
            Nj_minus_i = min_abundance_matrix.loc[i, j]
            Nj_star = max_abundance[j]

            if not np.isnan(cij) and not np.isnan(nij):
                numerator_ni += cij * Nj_minus_i * nij
                denominator_ni += cij * Nj_minus_i

            if not np.isnan(fij):
                fi_sum += fij * (Nj_minus_i / Nj_star)

        # Calculate Ni (Equation 10)
        if denominator_ni > 0:
            Ni = numerator_ni / denominator_ni
        else:
            Ni = np.nan

        # Calculate Fi (Equation 11)
        Fi = fi_sum

        # Set negative Ni and Fi to NaN
        if Ni is not None and Ni < 0:
            Ni = np.nan
        if Fi is not None and Fi < 0:
            Fi = np.nan

        results.loc[i, 'Ni'] = Ni
        results.loc[i, 'Fi'] = Fi

        # Coexistence judgment (Equation 13)
        if not np.isnan(Ni) and not np.isnan(Fi) and Ni < 1:
            results.loc[i, 'Coexistence'] = Fi <= 1 / (1 - Ni)
        else:
            results.loc[i, 'Coexistence'] = np.nan  # Set to NaN rather than False

    return results


# Main program
def main():
    # Read species abundance data
    abundance_df = pd.read_excel('Abundscale.xlsx', index_col=0)

    # Read global interaction matrix
    global_beta = pd.read_excel('InterMatrix.xlsx', index_col=0)

    # List of environmental gradient files
    env_files = [
        'EnvMatrix_Depth.xlsx',
        'EnvMatrix_NH4.N.xlsx',
        'EnvMatrix_NO2.N.xlsx',
        'EnvMatrix_NO3.N.xlsx',
        'EnvMatrix_Oxygen.xlsx',
        'EnvMatrix_pH.xlsx',
        'EnvMatrix_PO4.P.xlsx',
        'EnvMatrix_Salinity.xlsx',
        'EnvMatrix_SiO3.Si.xlsx',
        'EnvMatrix_Temperature.xlsx',
        'EnvMatrix_Turbidity.xlsx'
    ]

    # Calculate Ni, Fi, and coexistence judgment for each environmental gradient
    env_results = {}
    for env_file in env_files:
        print(f"Processing {env_file}...")
        env_beta = pd.read_excel(env_file, index_col=0)
        env_name = os.path.splitext(env_file)[0].replace('EnvMatrix_', '')
        env_results[env_name] = calculate_coexistence_metrics(env_beta, abundance_df)

        # Save results to Excel
        output_file = f"Coexistence_Results_{env_name}.xlsx"
        env_results[env_name].to_excel(output_file)
        print(f"Saved: {output_file}")

    # Calculate results for the global interaction matrix
    print("Processing global interaction matrix...")
    global_results = calculate_coexistence_metrics(global_beta, abundance_df)
    global_results.to_excel("Coexistence_Results_Global.xlsx")
    print("Saved: Coexistence_Results_Global.xlsx")

    # Calculate coexistence probability (Equation 14)
    coexistence_prob = pd.DataFrame(index=abundance_df.columns, columns=['Coexistence_Probability'])

    for species in abundance_df.columns:
        count_coexistence = 0
        total_valid = 0  # Number of valid points (Ni and Fi both non‑negative and non‑NaN)

        for env_name, result_df in env_results.items():
            if species in result_df.index:
                Ni = result_df.loc[species, 'Ni']
                Fi = result_df.loc[species, 'Fi']
                coexistence = result_df.loc[species, 'Coexistence']

                # Consider only points where Ni and Fi are both non‑negative and non‑NaN
                if not pd.isna(Ni) and not pd.isna(Fi) and Ni >= 0 and Fi >= 0:
                    total_valid += 1
                    if not pd.isna(coexistence) and coexistence:
                        count_coexistence += 1

        if total_valid > 0:
            coexistence_prob.loc[species, 'Coexistence_Probability'] = count_coexistence / total_valid
        else:
            coexistence_prob.loc[species, 'Coexistence_Probability'] = np.nan

    # Save coexistence probability results
    coexistence_prob.to_excel("Coexistence_Probability.xlsx")
    print("Saved: Coexistence_Probability.xlsx")

    print("All calculations completed!")


if __name__ == "__main__":
    main()