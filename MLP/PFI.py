"""
Task Type: Correlation Analysis
Usage Method: Deep Learning Models in Machine Learning
"""

'''
This code implements reading data from Excel, training a multi-layer perceptron (MLP) neural network model,
and determining the impact of each feature on the target variable through Permutation Feature Importance (PFI) analysis.
It also calculates model performance metrics: R², RMSE (original units), and NRMSE (normalized RMSE),
and outputs the results as Excel files and charts.
'''

import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.preprocessing import StandardScaler
import os

# ==================== Configuration Area (Modify according to actual situation) ====================
input_excel_path = "ALL_3_seasons.xlsx"          # Input Excel file path
# Other variables, 'Thalassiosira', 'Chaetoceros', 'Gymnodinium', 'Karenia mikimotoi', 'Chaetoceros tortissimus',
#                'Neoceratium fusus', 'Prorocentrum minimum', 'Chaetoceros curvisetus', 'Akashiwo sanguinea', 'Alexandrium',
#                'Gyrodinium dominans', 'Protoperidinium', 'Skeletonema costatum', 'ALL.cells'
target_vars = ['HAB']                             # Target variable(s) (multiple allowed)
features = [
    'Depth', 'NO3.N', 'PO4.P', 'NO2.N', 'NH4.N',
    'SiO3.Si', 'Temperature', 'Salinity', 'Oxygen',
    'Turbidity', 'pH'
]

# Feature labels (for plotting)
feature_labels = {
    'Depth': 'Depth',
    'NO3.N': r'NO$_3^-$-N',
    'PO4.P': r'PO$_4^{3-}$-P',
    'NO2.N': r'NO$_2^-$-N',
    'NH4.N': r'NH$_4^+$-N',
    'SiO3.Si': r'SiO$_3^{2-}$-Si',
    'Temperature': 'Temperature',
    'Salinity': 'Salinity',
    'Oxygen': 'DO',
    'Turbidity': 'Turbidity',
    'pH': 'pH'
}

# Feature colors (pastel color scheme, for plotting)
feature_colors = {
    'Depth': '#9BB7D4',
    'NO3.N': '#C3D69B',
    'PO4.P': '#F7D08A',
    'NO2.N': '#F4989C',
    'NH4.N': '#B2A1C7',
    'SiO3.Si': '#F9A875',
    'Temperature': '#F5B0CB',
    'Salinity': '#D7C49E',
    'Oxygen': '#87D7EB',
    'Turbidity': '#C4C4C4',
    'pH': '#F08080'
}

# Output path
base_output_dir = r"D:\AA\A小论文材料\PythonProjectMLP"
plt.rcParams['font.family'] = 'Arial'

# Random seed
def set_seed(seed=42):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

set_seed(42)

# ==================== Define MLP Model ====================
class SimpleMLP(nn.Module):
    def __init__(self, input_dim):
        super(SimpleMLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, 64)
        self.fc2 = nn.Linear(64, 32)
        self.fc3 = nn.Linear(32, 1)

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        return self.fc3(x)

# ==================== Main Processing Workflow ====================
excel_data = pd.ExcelFile(input_excel_path)

for target_var in target_vars:
    print(f"Processing dependent variable: {target_var}")
    target_dir = os.path.join(base_output_dir, target_var.replace('/', '_'))
    os.makedirs(target_dir, exist_ok=True)
    output_excel = os.path.join(target_dir, 'results.xlsx')
    output_figure = os.path.join(target_dir, 'importance.png')
    output_figure_pdf = os.path.join(target_dir, 'importance.pdf')

    results_dict = {}          # Store feature importance for each sheet
    performance_dict = {}       # Store performance metrics for each sheet

    for sheet_name in excel_data.sheet_names:
        set_seed(42)            # Fix seed for each sheet

        # Read data
        data = pd.read_excel(input_excel_path, sheet_name=sheet_name)
        # Clean: replace illegal characters, convert to numeric, forward/backward fill
        data = data.replace(to_replace=['/', ' ', '-', '--', '#DIV/0!'], value=np.nan)
        data = data.apply(pd.to_numeric, errors='coerce')
        data = data.ffill().bfill()

        if data.empty or not set(features).issubset(data.columns) or target_var not in data.columns:
            print(f"Skipping sheet: {sheet_name}")
            continue

        # Prepare feature X and target y
        X = data[features].apply(pd.to_numeric, errors='coerce').ffill()
        y = pd.to_numeric(data[target_var], errors='coerce').ffill()

        # Normalize target variable to [0,1] (for training)
        y_min, y_max = y.min(), y.max()
        if y_max != y_min:
            y_scaled = (y - y_min) / (y_max - y_min)
        else:
            y_scaled = y * 0   # Constant column

        # Standardize features
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        X_tensor = torch.tensor(X_scaled, dtype=torch.float32)
        y_tensor = torch.tensor(y_scaled.values, dtype=torch.float32).view(-1, 1)

        # Build model
        model = SimpleMLP(X_tensor.shape[1])
        criterion = nn.MSELoss()
        optimizer = optim.Adam(model.parameters(), lr=0.001)

        # Train
        for epoch in range(500):
            model.train()
            optimizer.zero_grad()
            loss = criterion(model(X_tensor), y_tensor)
            loss.backward()
            optimizer.step()
            if (epoch + 1) % 50 == 0 or epoch == 0:
                print(f"[{target_var} - {sheet_name}] Epoch {epoch+1}/500, Loss: {loss.item():.6f}")

        model.eval()
        with torch.no_grad():
            baseline_loss = criterion(model(X_tensor), y_tensor).item()

        # ========== New: Calculate model performance metrics ==========
        with torch.no_grad():
            y_pred_scaled = model(X_tensor)
            # Reverse normalization to original scale
            y_pred_original = y_pred_scaled * (y_max - y_min) + y_min
            y_true_original = y_tensor * (y_max - y_min) + y_min

            # RMSE (original units)
            mse = torch.mean((y_true_original - y_pred_original) ** 2)
            rmse = torch.sqrt(mse)

            # R²
            ss_res = torch.sum((y_true_original - y_pred_original) ** 2)
            ss_tot = torch.sum((y_true_original - torch.mean(y_true_original)) ** 2)
            if ss_tot == 0:
                r2 = 1.0
            else:
                r2 = 1 - ss_res / ss_tot

            # NRMSE (normalized RMSE = RMSE / range)
            y_range = y_max - y_min
            if y_range != 0:
                nrmse = rmse / y_range
            else:
                nrmse = float('nan')

            performance_dict[sheet_name] = {
                'R2': r2.item(),
                'RMSE': rmse.item(),
                'NRMSE': nrmse.item() if not np.isnan(nrmse) else nrmse
            }
            print(f"[{target_var} - {sheet_name}] Performance - R²: {r2.item():.4f}, RMSE: {rmse.item():.6f}, NRMSE: {nrmse.item():.4f}")
        # ========== End performance calculation ==========

        # ========== PFI Feature Importance Analysis (100 permutations) ==========
        num_permutations = 100
        feature_importances = np.zeros(X_tensor.shape[1])

        for i in range(X_tensor.shape[1]):
            permuted_losses = []
            for p in range(num_permutations):
                torch.manual_seed(42 + i * 100 + p)  # Ensure different seed per permutation but overall reproducible
                X_permuted = X_tensor.clone()
                X_permuted[:, i] = X_tensor[torch.randperm(X_tensor.size(0)), i]
                with torch.no_grad():
                    loss = criterion(model(X_permuted), y_tensor).item()
                    permuted_losses.append(loss)

            # Use the average loss increase over 100 permutations as final importance
            avg_permuted_loss = np.mean(permuted_losses)
            feature_importances[i] = avg_permuted_loss - baseline_loss

        total = feature_importances.sum()
        relative = feature_importances / total if total != 0 else feature_importances
        results = pd.DataFrame({'Feature': features, 'Importance': relative})
        results_dict[sheet_name] = results

    # ========== Save results to Excel ==========
    with pd.ExcelWriter(output_excel) as writer:
        # Write feature importance for each sheet
        for name, df in results_dict.items():
            df.to_excel(writer, sheet_name=name, index=False)

        # Write model performance metrics
        if performance_dict:
            perf_df = pd.DataFrame.from_dict(performance_dict, orient='index')
            perf_df.index.name = 'Season'
            perf_df = perf_df[['R2', 'RMSE', 'NRMSE']]  # Ensure column order
            perf_df.to_excel(writer, sheet_name='Model_Performance')

    # ========== Plot: Feature importance bar chart ==========
    num_sheets = len(results_dict)
    cols = 2
    rows = (num_sheets + cols - 1) // cols
    fig = plt.figure(figsize=(15, 6 * rows))
    gs = GridSpec(rows, cols, figure=fig)

    for idx, (sheet_name, df) in enumerate(results_dict.items()):
        df_sorted = df.sort_values(by='Importance', ascending=False).reset_index(drop=True)
        f_sorted = df_sorted['Feature']
        i_sorted = df_sorted['Importance']
        colors = [feature_colors.get(f, '#B0B0B0') for f in f_sorted]

        row, col = divmod(idx, cols)
        ax = fig.add_subplot(gs[row, col])
        ax.barh(f_sorted, i_sorted, color=colors)
        ax.invert_yaxis()
        ax.set_yticks(range(len(f_sorted)))
        ax.set_yticklabels([feature_labels.get(f, f) for f in f_sorted], fontsize=16)
        ax.tick_params(axis='x', labelsize=16)
        ax.set_title(sheet_name.replace(',', ''), fontsize=18)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_figure, dpi=300)
    plt.savefig(output_figure_pdf, dpi=300, format='pdf', bbox_inches='tight')
    plt.show()
    print(f"Completed saving for dependent variable {target_var}.")