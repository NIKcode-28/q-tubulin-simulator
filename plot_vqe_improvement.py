
import csv
import matplotlib.pyplot as plt

structures = ["1JFF", "9WDA", "9WDB"]

plt.figure(figsize=(10, 6))

for structure in structures:

    filename = f"{structure}_convergence.csv"

    with open(filename, "r") as file:
        data = list(csv.DictReader(file))

    energies = [
        float(row["total_energy_Ha"])
        for row in data
    ]

    initial_energy = energies[0]
    best_energy = initial_energy

    improvements = []

    for energy in energies:

        # Record the lowest energy found so far
        best_energy = min(best_energy, energy)

        # Improvement relative to the starting energy
        improvement = (initial_energy - best_energy) * 1000

        improvements.append(improvement)

    evaluations = range(1, len(energies) + 1)

    plt.plot(
        evaluations,
        improvements,
        marker="o",
        markersize=3,
        label=structure
    )

plt.xlabel("VQE Energy Evaluation")
plt.ylabel("Energy Improvement (mHa)")
plt.title("Energy Improvement During VQE Optimisation")

plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("vqe_energy_improvement.png", dpi=300)
plt.show()
