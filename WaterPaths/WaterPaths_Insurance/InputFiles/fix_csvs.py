import csv, os

files = [
    "rof_freq_table_state0_no_infra.csv",
    "rof_freq_table_state1_descoberto_infra.csv",
    "rof_freq_table_state2_tortoSM_infra.csv",
    "rof_freq_table_state3_both_infra.csv"
]

for fname in files:
    with open(fname, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = list(reader)

    thresholds = [row[0] for row in rows]
    freq_blocks = [row[1:] for row in rows]

    reversed_freqs = list(reversed(freq_blocks))

    corrected = []
    for i, thr in enumerate(thresholds):
        corrected.append([thr] + reversed_freqs[i])

    with open(fname, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(corrected)

    print(f"Fixed: {fname}")

print("\n--- state0 after fix (first 5 + last 5 rows) ---")
with open("rof_freq_table_state0_no_infra.csv", newline="") as f:
    reader = csv.reader(f)
    rows = list(reader)
for row in rows[:6]:   print(",".join(row))
print("...")
for row in rows[-5:]:
    print(",".join(row))
