import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── data ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(
    r"C:\Users\cmb247\repos\Abaqus\DCM_Scripting\BRC_graph.csv",
    sep="\t",
    header=0,
    names=["patient", "mJOA", "n_sites", "elem_id",
           "preload", "flexion", "extension", "comments"],
    usecols=["patient", "mJOA", "preload", "flexion", "extension"],
)
df = df.dropna(subset=["preload", "flexion", "extension"]).reset_index(drop=True)
df[["preload", "flexion", "extension"]] = df[["preload", "flexion", "extension"]].astype(float)
df["mJOA"] = df["mJOA"].astype(int)

patient_labels = ["P{:d}".format(int(p)) for p in df["patient"]]

conditions = [
    ("Preload",         "preload",    "o",  "#4C72B0"),
    ("Preload + Flexion",    "flexion",    "s",  "#DD8452"),
    ("Preload + Extension",  "extension",  "^",  "#55A868"),
]

# ── figure ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "Arial",
    "font.size":        11,
    "axes.linewidth":   0.8,
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "axes.spines.top":  False,
    "axes.spines.right": False,
})

fig, ax = plt.subplots(figsize=(5.5, 4.5))

# vertical line per patient connecting min to max strain across conditions
for _, row in df.iterrows():
    vals = [row["preload"], row["flexion"], row["extension"]]
    ax.vlines(row["mJOA"], min(vals), max(vals),
              color="#aaaaaa", linewidth=1.2, zorder=2)

for label, col, marker, colour in conditions:
    ax.scatter(df["mJOA"], df[col], label=label, color=colour,
               marker=marker, s=70, zorder=3)

# no-preload points for patient 1
p1_mjoa = int(df.loc[df["patient"] == 1, "mJOA"].iloc[0])
flexion_no_preload   = 0.04105915
extension_no_preload = 4.99996e-05

ax.vlines(p1_mjoa, extension_no_preload, flexion_no_preload,
          color="#aaaaaa", linewidth=1.2, zorder=2)
ax.scatter(p1_mjoa, flexion_no_preload, label="Flexion Only (No Preload)",
           marker="s", s=70, facecolors="none", edgecolors="#DD8452",
           linewidths=1.5, zorder=3)
ax.scatter(p1_mjoa, extension_no_preload, label="Extension Only (No Preload)",
           marker="^", s=70, facecolors="none", edgecolors="#55A868",
           linewidths=1.5, zorder=3)

# patient labels offset above each cluster of points
for i, (mjoa, row) in enumerate(df.iterrows()):
    y_max = df.loc[i, ["preload", "flexion", "extension"]].max()
    ax.annotate(patient_labels[i],
                xy=(df.loc[i, "mJOA"], y_max),
                xytext=(0, 7), textcoords="offset points",
                ha="center", fontsize=9, color="#444444")

ax.set_xlabel("mJOA Score")
ax.set_ylabel("Maximum Logarithmic Strain (LE)")
ax.set_ylim(-0.01, 0.35)
ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
#ax.grid(axis="y", linewidth=0.3, color="#cccccc", zorder=0)

# invert x so left = worse neurological function
#ax.invert_xaxis()
ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1] + 0.5)

ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.set_title("Maximum Single Site Strain Values Under 3 Load Types")
fig.tight_layout()
# fig.savefig(r"C:\Users\cmb247\repos\Abaqus\DCM_Scripting\BRC_graph.pdf", dpi=300)
fig.savefig(r"C:\Users\cmb247\repos\Abaqus\DCM_Scripting\BRC_graph_w_title_w+wo_preload.pdf", dpi=300)
# fig.savefig(r"C:\Users\cmb247\repos\Abaqus\DCM_Scripting\BRC_graph.png", dpi=300)

#plt.show()
print("Saved BRC_graph.png")
