# Fluid Flow Calculator

An interactive, competition-ready engineering decision-support tool for evaluating fluid flow hydraulics, friction factors, head losses, and pressure drops in circular pipe systems. Designed specifically for Petroleum, Mechanical, and Chemical Engineering applications, this tool applies rigorous fluid mechanics principles, robust numerical solvers for implicit equations, and clear unit management.

---

## Key Features

* **Comprehensive Fluid & Pipe Presets**: Built-in library of representative fluid properties (Water, Light Crude, Heavy Crude, Diesel, Gasoline, Ethylene Glycol) and standard Nominal Pipe Size (NPS) Schedule 40 dimensions, alongside full manual overrides.
* **Rigorous Friction Factor Calculation**:
  * **Laminar Flow ($Re < 2300$)**: Hagen-Poiseuille analytical equation ($f = 64/Re$).
  * **Transitional Flow ($2300 \le Re < 4000$)**: Calculates friction factor with explicit engineering caution regarding flow regime instability and empirical uncertainty.
  * **Turbulent Flow ($Re \ge 4000$)**: Colebrook-White implicit equation solved via Newton-Raphson iteration with Swamee-Jain explicit initial guess and numerical fallback safeguards.
* **Distinct Pressure Loss Breakdown**: Separates frictional wall losses ($\Delta P_f$) from hydrostatic elevation changes ($\Delta P_z$) with sign convention clarity ($\Delta z > 0$ for uphill incline, $\Delta z < 0$ for downhill gravity assistance).
* **Hydraulic Power Requirement**: Calculates hydraulic pumping power $P_{\text{hydraulic}} = \Delta P_{\text{total}} \cdot Q$ in kW and HP (explicitly distinguished from pump shaft power).
* **Interactive Engineering Visualizations**:
  * **Operating Sensitivity Curve**: Dynamic plot of frictional and total pressure loss versus volumetric flow rate $Q$.
  * **Interactive Moody Diagram**: Logarithmic $f$ vs. $Re$ grid with relative roughness ($\varepsilon/D$) iso-lines and current operating point coordinate.
* **Pandas Results Summary Table**: Real-time tabular output listing parameter symbols, calculated values, units, and engineering context with CSV export support.
* **Robust Input Guardrails**: Safeguards against non-physical inputs ($Q \le 0$, $D \le 0$, $L \le 0$, $\rho \le 0$, $\mu \le 0$, $\varepsilon < 0$) with warning badges for high flow velocities ($v > 10\text{ m/s}$).

---

## Governing Engineering Equations

### 1. Cross-Sectional Area ($A$)
$$A = \frac{\pi D^2}{4}$$
*where $D$ is the actual internal pipe diameter in meters.*

### 2. Average Flow Velocity ($v$)
$$v = \frac{Q}{A}$$
*where $Q$ is the volumetric flow rate in $\text{m}^3/\text{s}$.*

### 3. Reynolds Number ($Re$)
$$Re = \frac{\rho v D}{\mu}$$
*where $\rho$ is fluid density ($\text{kg/m}^3$) and $\mu$ is dynamic viscosity ($\text{Pa}\cdot\text{s}$).*

### 4. Darcy-Weisbach Friction Factor ($f$)
* **Laminar Flow ($Re < 2300$)**:
  $$f = \frac{64}{Re}$$
* **Turbulent Flow ($Re \ge 4000$)**: Implicit **Colebrook-White Equation**:
  $$\frac{1}{\sqrt{f}} = -2 \log_{10} \left( \frac{\varepsilon / D}{3.7} + \frac{2.51}{Re \sqrt{f}} \right)$$
  Solved numerically via Newton-Raphson method using **Swamee-Jain Explicit Equation** as initial guess:
  $$f_{\text{initial}} = \frac{0.25}{\left[ \log_{10} \left( \frac{\varepsilon/D}{3.7} + \frac{5.74}{Re^{0.9}} \right) \right]^2}$$
  *Convergence threshold: $|f_{k+1} - f_k| < 10^{-7}$ with maximum 50 iteration protection and fallback to Swamee-Jain.*
* **Transitional Flow ($2300 \le Re < 4000$)**:
  Friction factor evaluated using Colebrook-White / Swamee-Jain formulations with a prominent engineering caution flag highlighting physical flow regime instability and friction prediction uncertainty in the transitional zone.

### 5. Darcy-Weisbach Frictional Head Loss ($h_f$)
$$h_f = f \cdot \left(\frac{L}{D}\right) \cdot \left(\frac{v^2}{2g}\right)$$
*where $L$ is pipe length ($\text{m}$) and $g = 9.81\text{ m/s}^2$.*

### 6. Frictional Pressure Drop ($\Delta P_f$)
$$\Delta P_f = \rho g h_f = f \cdot \left(\frac{L}{D}\right) \cdot \left(\frac{\rho v^2}{2}\right)$$

### 7. Hydrostatic Elevation Pressure Change ($\Delta P_z$)
$$\Delta P_z = \rho g \Delta z$$
*Sign convention: $\Delta z = z_{\text{outlet}} - z_{\text{inlet}}$. Upward flow ($\Delta z > 0$) requires additional pressure; downhill flow ($\Delta z < 0$) gains hydrostatic head.*

### 8. Total System Pressure Drop Requirement ($\Delta P_{\text{total}}$)
$$\Delta P_{\text{total}} = \Delta P_f + \Delta P_z = \rho g (h_f + \Delta z)$$

### 9. Hydraulic Power Requirement ($P_{\text{hydraulic}}$)
$$P_{\text{hydraulic}} = \Delta P_{\text{total}} \cdot Q = \rho g Q (h_f + \Delta z)$$
*(Reported in Watts, kW, and Hydraulic Horsepower HP).*

---

## Inputs and Displayed Units

| Input Parameter | Internal SI Unit | Supported UI Units |
| :--- | :---: | :--- |
| Volumetric Flow Rate ($Q$) | $\text{m}^3/\text{s}$ | $\text{m}^3/\text{s}$, $\text{L/s}$, $\text{m}^3/\text{h}$, $\text{GPM}$ (gal/min), $\text{bbl/day}$ |
| Internal Pipe Diameter ($D$) | $\text{m}$ | $\text{mm}$, $\text{m}$, $\text{inches}$ |
| Pipe Length ($L$) | $\text{m}$ | $\text{m}$, $\text{km}$, $\text{feet}$ |
| Absolute Roughness ($\varepsilon$)| $\text{m}$ | $\text{mm}$, $\mu\text{m}$, $\text{inches}$ |
| Elevation Change ($\Delta z$) | $\text{m}$ | $\text{m}$, $\text{feet}$ |
| Fluid Density ($\rho$) | $\text{kg/m}^3$ | $\text{kg/m}^3$, $\text{g/cm}^3$, $\text{lb/ft}^3$, $^\circ\text{API}$ |
| Dynamic Viscosity ($\mu$) | $\text{Pa}\cdot\text{s}$ | $\text{Pa}\cdot\text{s}$, $\text{cP}$ (centipoise), $\text{mPa}\cdot\text{s}$ |

---

## Validation Rules & Safety Guards

* **Strict Input Validation**: If any critical parameter ($Q$, $D$, $L$, $\rho$, $\mu$) is $\le 0$ or roughness $\varepsilon < 0$, the application displays an explicit error message and halts calculation without crashing.
* **Velocity Alert ($v > 10\text{ m/s}$)**: Triggers an engineering warning highlighting potential risks of pipe wall erosion, acoustic noise, or pressure surging depending on pipe metallurgy and liquid phase.
* **Transitional Regime Alert ($2300 \le Re < 4000$)**: Highlights flow instability in the critical region.

---

## Technologies Used

* **Python 3.10+**
* **Streamlit**: Web interface and interactive layout framework
* **Pandas**: Tabular data structuring and CSV export functionality
* **NumPy**: Numerical computations and vectorization
* **Plotly**: Dynamic interactive charts (operating curve & Moody diagram)

---

## How to Run Locally

1. **Clone the repository**:
   ```bash
   git clone <GITHUB_URL>
   cd fluid-flow-calculator
   ```

2. **Set up a virtual environment (optional but recommended)**:
   ```bash
   python3 -m venv venv
   source venv/bin/activate   # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

4. **Launch the Streamlit app**:
   ```bash
   streamlit run app.py
   ```

5. **Open browser**: Navigate to `http://localhost:8501`.

---

## Deployment

This application is configured for seamless deployment on **Streamlit Community Cloud**:
1. Push source code to GitHub.
2. Connect repository to Streamlit Cloud dashboard.
3. Set main file path to `app.py`.

* **LIVE APP URL**: [TO BE ADDED AFTER DEPLOYMENT]
* **GITHUB REPOSITORY URL**: [TO BE ADDED AFTER DEPLOYMENT]
