"""
02_visualize.py
Produce one clear bar chart showing which variables most influence tip amount,
using standardized regression coefficients (so magnitudes are comparable
across variables measured in different units/scales).
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

df = pd.read_csv("tips_clean.csv")

# One-hot encode categoricals (drop_first avoids dummy trap)
X = pd.get_dummies(
    df[["total_bill", "size", "sex", "smoker", "day", "time"]],
    drop_first=True,
)
y = df["tip"]

# Standardize features so coefficients are directly comparable
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

model = LinearRegression().fit(X_scaled, y)

coef_df = pd.DataFrame({
    "feature": X.columns,
    "std_coef": model.coef_
}).sort_values("std_coef", key=abs, ascending=True)

# --- Bar chart ---
colors = ["#d62728" if c < 0 else "#2ca02c" for c in coef_df["std_coef"]]

fig, ax = plt.subplots(figsize=(8, 5))
ax.barh(coef_df["feature"], coef_df["std_coef"], color=colors)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Standardized coefficient (effect on tip, in $, per 1 SD change)")
ax.set_title("What Drives Tip Amount? Standardized Regression Coefficients")
ax.grid(axis="x", linestyle="--", alpha=0.4)

for i, v in enumerate(coef_df["std_coef"]):
    ax.text(v + (0.02 if v >= 0 else -0.02), i, f"{v:.2f}",
            va="center", ha="left" if v >= 0 else "right", fontsize=9)

plt.tight_layout()
plt.savefig("tip_drivers_barchart.png", dpi=150)
print("Saved chart -> tip_drivers_barchart.png")
print("\nStandardized coefficients (sorted by magnitude):")
print(coef_df.sort_values("std_coef", key=abs, ascending=False).to_string(index=False))
