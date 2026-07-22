"""
01_clean.py
Load and clean the tips dataset.
"""
import pandas as pd

df = pd.read_csv("tips.csv")

print("Shape:", df.shape)
print("\nDtypes:\n", df.dtypes)
print("\nMissing values:\n", df.isna().sum())
print("\nDuplicate rows:", df.duplicated().sum())
print("\nSummary stats:\n", df.describe(include="all"))

# --- Cleaning steps ---
# 1. Drop exact duplicate rows (none expected, but safe practice)
df = df.drop_duplicates()

# 2. Standardise categorical text (strip whitespace, consistent case)
cat_cols = ["sex", "smoker", "day", "time"]
for c in cat_cols:
    df[c] = df[c].astype(str).str.strip()

# 3. Sanity-check ranges — no negative bills/tips, size >= 1
df = df[(df["total_bill"] > 0) & (df["tip"] >= 0) & (df["size"] >= 1)]

# 4. Feature engineering: tip percentage (useful for the "influence" story)
df["tip_pct"] = df["tip"] / df["total_bill"] * 100

df.to_csv("tips_clean.csv", index=False)
print("\nCleaned data saved -> tips_clean.csv")
print("Final shape:", df.shape)
