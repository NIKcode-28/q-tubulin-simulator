
import csv
import matplotlib.pyplot as plt

# Classical reference energies (Hartree)
references = {
    "1JFF": -700.90120455,
    "9WDA": -700.96681343,
    "9WDB": -700.96121968
}

plt.figure(figsize=(10, 6))

for structure, reference in references.items():

    # Read the convergence CSV
    filename = f"{structure}_convergence.csv"

    with open(filename, "r") as file:
        reader = csv.DictReader(file)
        data = list(reader)

    evaluations = []
    errors = []

    best_energy = float("inf")

    for row in data:
        evaluation = int(row["evaluation"])
        energy = float(row["total_energy_Ha"])

        # Keep track of the lowest energy found so far
        best_energy = min(best_energy, energy)

        # Calculate error in millihartree
        error = abs(best_energy - reference) * 1000

        evaluations.append(evaluation)
        errors.append(error)

    # Plot the error for this structure
    plt.plot(
        evaluations,
        errors,
        marker="o",
        markersize=3,
        label=structure
    )

# Chemical accuracy reference
plt.axhline(
    y=1.6,
    color="red",
    linestyle="--",
    label="Chemical accuracy (1.6 mHa)"
)

plt.xlabel("VQE Energy Evaluation")
plt.ylabel("Absolute Energy Error (mHa)")
plt.title("VQE Convergence for Tubulin Molecular Models")
plt.yscale("log")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("vqe_convergence_comparison.png", dpi=300)
plt.show()
