# SSMK-MPC: Data-Driven Micromotion Control for Soft Robotic Laser Surgery  

This repository provides the implementation of an **MPC controller with a data-driven surrogate model** for dynamic micromotion control in laser surgery.  

---

## ✨ Key Features

- **Surrogate Modeling**  
  Derived by adapting **Extended Dynamic Mode Decomposition (EDMD)** to approximate the Koopman operator on **low-dimensional Spectral Submanifold (SSM) coordinates**.  

- **Efficiency & Accuracy**  
  Captures task-relevant dynamics with high fidelity while preserving linearity for efficient **MPC implementation at 100 Hz**.  

- **Experimental Validation**  
  Demonstrated on **millimeter-scale figure-8 trajectory tracking**, where **SSMK-MPC** outperforms both SSMP and Koopman models, effectively suppressing oscillations seen in direct Koopman-based approaches.  

- **Surgical Relevance**  
  Precise micromotion control enables sharp incisions, reduced collateral thermal damage, and motion compensation for respiration and heartbeat.  

---

## 🚀 Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/your-username/ssmk-mpc.git
cd ssmk-mpc

<img width="1666" height="835" alt="Architecture" src="https://github.com/user-attachments/assets/1e8cb86f-07f9-498d-a345-6e32f0343d7c" />
