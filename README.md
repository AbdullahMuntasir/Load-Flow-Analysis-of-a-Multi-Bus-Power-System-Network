# Load Flow Analysis of a 9-Bus Power System Network

 📌 Overview
This project implements and compares **Load Flow (Power Flow) Analysis** of a **9-bus power system** using:
- **MATLAB** (Gauss–Seidel Method)
- **PowerWorld Simulator** (Newton–Raphson Method)

The aim is to determine **bus voltages, voltage angles, line power flows, and system losses**, and to analyze differences in results between the two methods.

---

📖 Background
Load Flow Analysis is a critical tool in power system engineering for:
- Determining voltage magnitude and phase angle at each bus
- Calculating active (MW) and reactive (MVAr) power flows
- Evaluating system losses and voltage regulation
- Ensuring stable and efficient grid operation

This project models a **9-bus test system** with:
- 1 Slack Bus
- 2 PV (Generator) Buses
- 6 PQ (Load) Buses
- 10 Transmission Lines

---

 🛠 Tools & Methods
 1. MATLAB**
- Method: **Gauss–Seidel Iterative Method**
- Custom implementation with per-unit conversion
- Includes PV bus reactive power limit checks
- Outputs:
  - Bus voltages & angles (p.u. & kV)
  - Line flows (MW, MVAr)
  - System losses
  - Graphs for voltage profile & line flows

 2. PowerWorld Simulator**
- Method: **Newton–Raphson**
- Graphical single-line diagram
- Built-in load flow solver
- Outputs:
  - Tabular bus voltage and flow data
  - Loss reports
  - Voltage profile visualization
