"""
===============================================================================
AI-ASSISTED DEVELOPMENT DOCUMENTATION
===============================================================================
AI Tools Used:
- Google AI Studio Build Environment (Gemini 3.6 Flash model)
- Python 3.10+ with Streamlit, Pandas, NumPy, and Plotly

Key Prompts Used During Development:
1. "Create a competition-ready Streamlit fluid flow calculator for petroleum and
   mechanical engineering applications with strict SI calculations, Colebrook-White
   friction factor numerical solver, and input validation."
2. "Implement robust Newton-Raphson iteration for the Colebrook-White equation with
   Swamee-Jain explicit initial guess, 50-iteration safeguard, fallback mechanisms,
   and distinct frictional versus elevation pressure components."
3. "Design interactive Plotly sensitivity charts and Moody diagram alongside a Pandas
   summary table, unit converters, and balanced engineering warning badges."

Manual Verification & Corrections Performed:
- Verified Colebrook-White Newton-Raphson solver against authoritative textbook Moody
  diagram values (e.g., Re = 10^5, relative roughness e/D = 0.001 yields f ≈ 0.0222).
- Verified Hagen-Poiseuille analytical solution (f = 64/Re) for laminar flow (Re = 1000
  yields f = 0.0640).
- Manually checked sign conventions for elevation changes: positive delta_z (uphill flow)
  increases required inlet pressure (delta_P_z > 0), while negative delta_z (downhill flow)
  assists flow (delta_P_z < 0).
- Verified unit conversion constants (1 GPM = 6.30902e-5 m³/s; 1 cP = 1e-3 Pa·s;
  1 inch = 0.0254 m; 1 bbl/day = 1.84013e-6 m³/s).
- Refined flow velocity warning wording to avoid declaring any single threshold as a
  universal erosion limit, noting that erosion depends on metallurgy, liquid phase,
  and industry standards.
===============================================================================
"""

import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# -----------------------------------------------------------------------------
# PAGE CONFIGURATION & STYLING
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="Fluid Flow Calculator",
    page_icon="⛽",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for polished engineering aesthetics
st.markdown(
    """
    <style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.05rem;
        color: #4B5563;
        margin-bottom: 1.2rem;
    }
    .metric-card {
        background-color: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 12px;
        text-align: center;
    }
    .metric-value {
        font-size: 1.4rem;
        font-weight: 700;
        color: #0F172A;
    }
    .metric-label {
        font-size: 0.85rem;
        color: #64748B;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# -----------------------------------------------------------------------------
# CONSTANTS & DATABASE PRESETS
# -----------------------------------------------------------------------------
G_ACCEL = 9.81  # Gravitational acceleration, m/s²

FLUID_PRESETS = {
    "Water (20°C)": {
        "rho": 998.2,  # kg/m³
        "mu": 0.001002,  # Pa·s
        "note": "Standard fresh water properties at atmospheric pressure.",
    },
    "Light Crude Oil (35°API)": {
        "rho": 850.0,
        "mu": 0.007000,
        "note": "Representative properties at reference temperature (15.6°C / 60°F).",
    },
    "Heavy Crude Oil (20°API)": {
        "rho": 934.0,
        "mu": 0.100000,
        "note": "Representative viscous crude properties at reference temperature.",
    },
    "Diesel Fuel": {
        "rho": 832.0,
        "mu": 0.003000,
        "note": "Representative commercial diesel fuel properties.",
    },
    "Gasoline": {
        "rho": 740.0,
        "mu": 0.000600,
        "note": "Representative motor gasoline fuel properties.",
    },
    "Ethylene Glycol (100% solution)": {
        "rho": 1113.0,
        "mu": 0.016000,
        "note": "Standard industrial ethylene glycol heat transfer fluid.",
    },
    "Custom Fluid": {
        "rho": 1000.0,
        "mu": 0.001000,
        "note": "User-defined fluid density and dynamic viscosity.",
    },
}

PIPE_SCHEDULE_PRESETS = {
    "Custom Diameter": None,
    'NPS 1" Sch 40 (ID: 26.64 mm / 1.049 in)': 26.64,
    'NPS 2" Sch 40 (ID: 52.50 mm / 2.067 in)': 52.50,
    'NPS 3" Sch 40 (ID: 77.92 mm / 3.068 in)': 77.92,
    'NPS 4" Sch 40 (ID: 102.26 mm / 4.026 in)': 102.26,
    'NPS 6" Sch 40 (ID: 154.08 mm / 6.065 in)': 154.08,
    'NPS 8" Sch 40 (ID: 202.72 mm / 7.981 in)': 202.72,
    'NPS 10" Sch 40 (ID: 254.46 mm / 10.020 in)': 254.46,
    'NPS 12" Sch 40 (ID: 303.22 mm / 11.938 in)': 303.22,
}

PIPE_ROUGHNESS_PRESETS = {
    "Commercial Steel / Wrought Iron (0.045 mm)": 0.045,
    "Smooth Plastic / PVC / Glass (0.0015 mm)": 0.0015,
    "Drawn Brass / Copper / Tubing (0.0025 mm)": 0.0025,
    "Stainless Steel (0.015 mm)": 0.015,
    "Galvanized Iron (0.15 mm)": 0.15,
    "Cast Iron (0.26 mm)": 0.26,
    "Concrete (1.0 mm)": 1.00,
    "Custom Roughness": None,
}

# -----------------------------------------------------------------------------
# HELPER CALCULATION FUNCTIONS
# -----------------------------------------------------------------------------


def calculate_area(diameter_m: float) -> float:
    """Calculates cross-sectional hydraulic area A = pi * D^2 / 4."""
    return math.pi * (diameter_m**2) / 4.0


def calculate_velocity(q_m3s: float, area_m2: float) -> float:
    """Calculates average flow velocity v = Q / A."""
    return q_m3s / area_m2


def calculate_reynolds_number(
    rho: float, velocity: float, diameter_m: float, mu: float
) -> float:
    """Calculates dimensionless Reynolds Number Re = (rho * v * D) / mu."""
    return (rho * velocity * diameter_m) / mu


def calculate_swamee_jain_friction_factor(
    reynolds: float, rel_roughness: float
) -> float:
    """
    Calculates explicit initial guess for friction factor using Swamee-Jain (1976):
    f = 0.25 / [log10( (e/D)/3.7 + 5.74 / (Re^0.9) )]^2
    """
    term = (rel_roughness / 3.7) + (5.74 / (reynolds**0.9))
    log_term = math.log10(term)
    return 0.25 / (log_term**2)


def calculate_colebrook_white_friction_factor(
    reynolds: float, rel_roughness: float, max_iter: int = 50, tol: float = 1e-7
) -> tuple[float, bool, int]:
    """
    Solves implicit Colebrook-White equation for turbulent flow via Newton-Raphson:
    g(f) = 1/sqrt(f) + 2*log10( (e/D)/3.7 + 2.51 / (Re*sqrt(f)) ) = 0

    Returns:
        (f_value, converged_flag, iterations_used)
    """
    # Use Swamee-Jain as explicit initial guess
    f_curr = calculate_swamee_jain_friction_factor(reynolds, rel_roughness)

    for i in range(1, max_iter + 1):
        sqrt_f = math.sqrt(f_curr)
        arg = (rel_roughness / 3.7) + (2.51 / (reynolds * sqrt_f))

        if arg <= 0 or f_curr <= 0:
            break

        # Residual g(f)
        g_val = (1.0 / sqrt_f) + 2.0 * math.log10(arg)

        # Derivative g'(f)
        dg_df = -0.5 * (f_curr ** (-1.5)) - (
            2.51 / (math.log(10.0) * reynolds * (f_curr**1.5) * arg)
        )

        f_next = f_curr - (g_val / dg_df)

        # Numerical safeguard against non-positive values
        if f_next <= 0:
            f_next = f_curr / 2.0

        if abs(f_next - f_curr) < tol:
            return f_next, True, i

        f_curr = f_next

    # Fallback to Swamee-Jain explicit approximation if Newton-Raphson exceeds max_iter
    f_fallback = calculate_swamee_jain_friction_factor(reynolds, rel_roughness)
    return f_fallback, False, max_iter


def calculate_friction_factor(
    reynolds: float, rel_roughness: float
) -> tuple[float, str, str]:
    """
    Determines friction factor f and regime based on Reynolds number:
    - Re < 2300: Laminar (f = 64 / Re)
    - 2300 <= Re < 4000: Critical Transitional Regime (Colebrook-White with warning)
    - Re >= 4000: Fully Turbulent (Colebrook-White Newton-Raphson)

    Returns:
        (f_value, regime_name, calculation_method_description)
    """
    if reynolds < 2300.0:
        f_val = 64.0 / reynolds
        return f_val, "Laminar", "Hagen-Poiseuille Analytical (f = 64/Re)"
    elif 2300.0 <= reynolds < 4000.0:
        f_val, converged, iters = calculate_colebrook_white_friction_factor(
            reynolds, rel_roughness
        )
        method_str = (
            f"Colebrook-White Solver ({iters} iter) [Critical Transitional Region]"
        )
        return f_val, "Transitional", method_str
    else:
        f_val, converged, iters = calculate_colebrook_white_friction_factor(
            reynolds, rel_roughness
        )
        method_str = f"Colebrook-White Newton-Raphson Solver ({iters} iter)"
        return f_val, "Turbulent", method_str


def calculate_head_loss(
    f: float, length_m: float, diameter_m: float, velocity: float
) -> float:
    """Calculates Darcy-Weisbach frictional head loss h_f = f * (L/D) * (v^2 / (2g))."""
    return f * (length_m / diameter_m) * ((velocity**2) / (2.0 * G_ACCEL))


def calculate_pressure_drops(
    rho: float, h_f: float, delta_z: float
) -> tuple[float, float, float]:
    """
    Calculates pressure drops:
    - Frictional drop: delta_P_f = rho * g * h_f (Pa)
    - Elevation change: delta_P_z = rho * g * delta_z (Pa)
    - Total required system drop: delta_P_total = delta_P_f + delta_P_z (Pa)
    """
    delta_P_f = rho * G_ACCEL * h_f
    delta_P_z = rho * G_ACCEL * delta_z
    delta_P_total = delta_P_f + delta_P_z
    return delta_P_f, delta_P_z, delta_P_total


def calculate_hydraulic_power(delta_P_total_pa: float, q_m3s: float) -> tuple[float, float, float]:
    """
    Calculates hydraulic power requirement P_hydraulic = delta_P_total * Q.
    Returns (Watts, kW, HP).
    """
    p_watts = delta_P_total_pa * q_m3s
    p_kw = p_watts / 1000.0
    p_hp = p_watts / 745.7  # 1 HP = 745.7 W
    return p_watts, p_kw, p_hp


# -----------------------------------------------------------------------------
# MAIN APPLICATION INTERFACE
# -----------------------------------------------------------------------------

st.markdown('<div class="main-header">FLUID FLOW CALCULATOR</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="sub-header">An interactive university engineering decision-support tool for evaluating fluid flow, friction factors, and pressure losses in circular pipes.</div>',
    unsafe_allow_html=True,
)

# User Instructions & Assumptions Expander
with st.expander("📌 User Instructions & Engineering Operating Guidelines", expanded=False):
    st.markdown(
        """
        1. **Select or Input Fluid Properties**: Choose from petroleum/industrial presets or select *Custom Fluid* to input density ($\rho$) and viscosity ($\mu$).
        2. **Specify Pipe System Parameters**: Select Nominal Pipe Size (NPS) Schedule 40 or enter actual internal diameter ($D$), pipe length ($L$), surface roughness ($\varepsilon$), and elevation change ($\Delta z$).
        3. **Set Volumetric Flow Rate**: Input flow rate ($Q$) in preferred units ($\text{m}^3/\text{s}$, $\text{L/s}$, $\text{m}^3/\text{h}$, $\text{GPM}$, or $\text{bbl/day}$).
        4. **Review Hydraulics & Metrics**: Examine flow velocity ($v$), Reynolds number ($Re$), friction factor ($f$), frictional head loss ($h_f$), and distinct pressure drops.
        5. **Analyze Sensitivity & Charts**: Explore interactive operating curves ($\Delta P$ vs. $Q$) and the interactive Moody diagram.
        """
    )

with st.expander("📚 Engineering Assumptions & Physical Limitations", expanded=False):
    st.markdown(
        """
        * **Steady-State Incompressible Flow**: Fluid density remains constant along pipe length; temperature and pressure effects on fluid properties are neglected.
        * **Single-Phase Newtonian Fluids**: Equations apply to homogeneous Newtonian fluids with constant viscosity.
        * **Circular Pipe Geometry**: Calculations assume clean, circular cross-sections with uniform internal diameter $D$ and uniform absolute wall roughness $\varepsilon$.
        * **Major Frictional Losses Focus**: Head loss calculations evaluate straight pipe major friction loss ($h_f$). Minor losses from valves/fittings are excluded in this baseline model.
        * **Elevation Sign Convention**: Positive elevation change ($\Delta z > 0$) represents uphill flow requiring additional pump pressure head. Negative elevation ($\Delta z < 0$) represents downhill flow gaining hydrostatic head.
        * **Petroleum Presets**: Preset properties for crude oils and petroleum fuels are representative values at reference conditions ($15.6^\circ\text{C} / 60^\circ\text{F}$) and do not replace laboratory assay measurements for specific fluid batches.
        """
    )

# -----------------------------------------------------------------------------
# SIDEBAR CONTROLS
# -----------------------------------------------------------------------------
st.sidebar.header("⚙️ System Input Parameters")

# --- 1. FLUID SELECTION ---
st.sidebar.subheader("1. Fluid Properties")
fluid_choice = st.sidebar.selectbox("Fluid Preset", list(FLUID_PRESETS.keys()), index=0)
fluid_info = FLUID_PRESETS[fluid_choice]
st.sidebar.caption(f"ℹ️ {fluid_info['note']}")

col_rho_unit, col_mu_unit = st.sidebar.columns(2)

# Density Units
rho_unit = col_rho_unit.selectbox("Density Unit", ["kg/m³", "g/cm³", "lb/ft³", "°API"], index=0)
default_rho = fluid_info["rho"]

if rho_unit == "g/cm³":
    default_rho = default_rho / 1000.0
elif rho_unit == "lb/ft³":
    default_rho = default_rho * 0.062428
elif rho_unit == "°API":
    # 141.5 / (SG) - 131.5 => SG = 141.5 / (°API + 131.5)
    sg = default_rho / 999.0
    default_rho = (141.5 / sg) - 131.5

input_rho_val = st.sidebar.number_input(
    f"Density ({rho_unit})", value=float(default_rho), format="%.4f"
)

# Convert density to kg/m³
if rho_unit == "kg/m³":
    rho_kg_m3 = input_rho_val
elif rho_unit == "g/cm³":
    rho_kg_m3 = input_rho_val * 1000.0
elif rho_unit == "lb/ft³":
    rho_kg_m3 = input_rho_val / 0.062428
elif rho_unit == "°API":
    if input_rho_val > -131.5:
        sg_calc = 141.5 / (input_rho_val + 131.5)
        rho_kg_m3 = sg_calc * 999.0
    else:
        rho_kg_m3 = 1000.0

# Viscosity Units
mu_unit = col_mu_unit.selectbox("Viscosity Unit", ["Pa·s", "mPa·s / cP"], index=0)
default_mu = fluid_info["mu"]
if mu_unit == "mPa·s / cP":
    default_mu = default_mu * 1000.0

input_mu_val = st.sidebar.number_input(
    f"Dynamic Viscosity ({mu_unit})", value=float(default_mu), format="%.6f"
)

# Convert viscosity to Pa·s
if mu_unit == "Pa·s":
    mu_pas = input_mu_val
else:
    mu_pas = input_mu_val / 1000.0


# --- 2. FLOW RATE SELECTION ---
st.sidebar.subheader("2. Flow Rate")
q_unit = st.sidebar.selectbox(
    "Flow Rate Unit", ["m³/s", "L/s", "m³/h", "GPM (gal/min)", "bbl/day"], index=1
)

# Default flow rate = 10 L/s = 0.01 m³/s
if q_unit == "m³/s":
    default_q = 0.010
elif q_unit == "L/s":
    default_q = 10.0
elif q_unit == "m³/h":
    default_q = 36.0
elif q_unit == "GPM (gal/min)":
    default_q = 158.5
elif q_unit == "bbl/day":
    default_q = 5434.4

input_q_val = st.sidebar.number_input(f"Volumetric Flow Rate Q ({q_unit})", value=float(default_q), format="%.4f")

# Convert Q to m³/s
if q_unit == "m³/s":
    q_m3s = input_q_val
elif q_unit == "L/s":
    q_m3s = input_q_val / 1000.0
elif q_unit == "m³/h":
    q_m3s = input_q_val / 3600.0
elif q_unit == "GPM (gal/min)":
    q_m3s = input_q_val * 6.30902e-5
elif q_unit == "bbl/day":
    q_m3s = input_q_val * 1.84013e-6


# --- 3. PIPE & SYSTEM GEOMETRY ---
st.sidebar.subheader("3. Pipe Geometry & Material")

pipe_preset_choice = st.sidebar.selectbox("Pipe Size Preset", list(PIPE_SCHEDULE_PRESETS.keys()), index=2)
preset_diameter_mm = PIPE_SCHEDULE_PRESETS[pipe_preset_choice]

dia_unit = st.sidebar.selectbox("Diameter Unit", ["mm", "m", "inches"], index=0)

if preset_diameter_mm is not None:
    if dia_unit == "mm":
        default_dia = preset_diameter_mm
    elif dia_unit == "m":
        default_dia = preset_diameter_mm / 1000.0
    elif dia_unit == "inches":
        default_dia = preset_diameter_mm / 25.4
else:
    default_dia = 52.50 if dia_unit == "mm" else (0.0525 if dia_unit == "m" else 2.067)

input_dia_val = st.sidebar.number_input(
    f"Actual Internal Diameter D ({dia_unit})", value=float(default_dia), format="%.4f"
)
st.sidebar.caption("⚠️ Note: Always specify actual *internal* diameter D for accurate hydraulic area calculations.")

# Convert Diameter to meters
if dia_unit == "mm":
    diameter_m = input_dia_val / 1000.0
elif dia_unit == "m":
    diameter_m = input_dia_val
elif dia_unit == "inches":
    diameter_m = input_dia_val * 0.0254

# Pipe Length
len_unit = st.sidebar.selectbox("Pipe Length Unit", ["m", "km", "feet"], index=0)
default_len = 100.0 if len_unit == "m" else (0.1 if len_unit == "km" else 328.08)
input_len_val = st.sidebar.number_input(f"Pipe Length L ({len_unit})", value=float(default_len), format="%.2f")

if len_unit == "m":
    length_m = input_len_val
elif len_unit == "km":
    length_m = input_len_val * 1000.0
elif len_unit == "feet":
    length_m = input_len_val * 0.3048

# Pipe Roughness
rough_preset_choice = st.sidebar.selectbox("Material Roughness Preset", list(PIPE_ROUGHNESS_PRESETS.keys()), index=0)
preset_rough_mm = PIPE_ROUGHNESS_PRESETS[rough_preset_choice]

rough_unit = st.sidebar.selectbox("Roughness Unit", ["mm", "µm", "inches"], index=0)

if preset_rough_mm is not None:
    if rough_unit == "mm":
        default_rough = preset_rough_mm
    elif rough_unit == "µm":
        default_rough = preset_rough_mm * 1000.0
    elif rough_unit == "inches":
        default_rough = preset_rough_mm / 25.4
else:
    default_rough = 0.045 if rough_unit == "mm" else (45.0 if rough_unit == "µm" else 0.00177)

input_rough_val = st.sidebar.number_input(
    f"Absolute Surface Roughness ε ({rough_unit})", value=float(default_rough), format="%.6f"
)

if rough_unit == "mm":
    roughness_m = input_rough_val / 1000.0
elif rough_unit == "µm":
    roughness_m = input_rough_val / 1e6
elif rough_unit == "inches":
    roughness_m = input_rough_val * 0.0254

# Elevation Change
elev_unit = st.sidebar.selectbox("Elevation Unit", ["m", "feet"], index=0)
input_elev_val = st.sidebar.number_input(
    f"Elevation Change Δz ({elev_unit})", value=0.0, format="%.2f"
)
st.sidebar.caption("Sign convention: Δz > 0 for uphill flow; Δz < 0 for downhill flow.")

if elev_unit == "m":
    delta_z_m = input_elev_val
else:
    delta_z_m = input_elev_val * 0.3048


# -----------------------------------------------------------------------------
# INPUT VALIDATION & GUARDRAILS
# -----------------------------------------------------------------------------
if q_m3s <= 0.0:
    st.error("⛔ **Invalid Input Error**: Volumetric flow rate Q must be strictly positive ($Q > 0$). Please correct the flow rate input.")
    st.stop()

if diameter_m <= 0.0:
    st.error("⛔ **Invalid Input Error**: Pipe internal diameter D must be strictly positive ($D > 0$). Please correct the diameter input.")
    st.stop()

if length_m <= 0.0:
    st.error("⛔ **Invalid Input Error**: Pipe length L must be strictly positive ($L > 0$). Please correct the length input.")
    st.stop()

if rho_kg_m3 <= 0.0:
    st.error("⛔ **Invalid Input Error**: Fluid density ρ must be strictly positive ($\rho > 0$). Please correct density input.")
    st.stop()

if mu_pas <= 0.0:
    st.error("⛔ **Invalid Input Error**: Dynamic viscosity μ must be strictly positive ($\mu > 0$). Please correct viscosity input.")
    st.stop()

if roughness_m < 0.0:
    st.error("⛔ **Invalid Input Error**: Absolute pipe roughness ε cannot be negative ($\varepsilon \ge 0$). Please correct roughness input.")
    st.stop()


# -----------------------------------------------------------------------------
# CORE ENGINEERING COMPUTATIONS
# -----------------------------------------------------------------------------
area_m2 = calculate_area(diameter_m)
velocity_m_s = calculate_velocity(q_m3s, area_m2)
reynolds = calculate_reynolds_number(rho_kg_m3, velocity_m_s, diameter_m, mu_pas)
relative_roughness = roughness_m / diameter_m

f_factor, flow_regime, method_desc = calculate_friction_factor(reynolds, relative_roughness)
head_loss_m = calculate_head_loss(f_factor, length_m, diameter_m, velocity_m_s)

delta_P_f_pa, delta_P_z_pa, delta_P_total_pa = calculate_pressure_drops(
    rho_kg_m3, head_loss_m, delta_z_m
)

p_watts, p_kw, p_hp = calculate_hydraulic_power(delta_P_total_pa, q_m3s)


# -----------------------------------------------------------------------------
# ENGINEERING CAUTION NOTIFICATIONS
# -----------------------------------------------------------------------------
if flow_regime == "Transitional":
    st.warning(
        f"⚠️ **Critical Transitional Flow Regime Alert ($2300 \le Re < 4000$)**: Calculated Reynolds number ($Re = {reynolds:,.0f}$) falls in the critical transitional zone. Flow exhibits temporal instability, and friction factors carry high empirical uncertainty."
    )

if velocity_m_s > 10.0:
    st.warning(
        f"⚠️ **High Velocity Alert ($v = {velocity_m_s:.2f} \\text{{ m/s}}$)**: Calculated flow velocity is relatively high ($v > 10 \\text{{ m/s}}$). Depending on pipe metallurgy, liquid phase, and industry standards (e.g. API RP 14E for petroleum piping), high velocities may increase risks of acoustic noise, hydraulic surge, or pipe wall erosion."
    )


# -----------------------------------------------------------------------------
# RESULTS DISPLAY — METRIC CARDS OVERVIEW
# -----------------------------------------------------------------------------
st.subheader("📊 Hydraulic Results Overview")

m_col1, m_col2, m_col3, m_col4 = st.columns(4)

with m_col1:
    st.metric(label="Flow Velocity (v)", value=f"{velocity_m_s:.2f} m/s", delta=f"Area: {area_m2*1e4:.2f} cm²")

with m_col2:
    st.metric(label="Reynolds Number (Re)", value=f"{reynolds:,.0f}", delta=f"Regime: {flow_regime}")

with m_col3:
    st.metric(label="Darcy Friction Factor (f)", value=f"{f_factor:.5f}", delta=f"ε/D = {relative_roughness:.5f}")

with m_col4:
    st.metric(label="Head Loss (h_f)", value=f"{head_loss_m:.2f} m", delta=f"Length: {length_m:.1f} m")

m_col5, m_col6, m_col7, m_col8 = st.columns(4)

with m_col5:
    st.metric(label="Frictional ΔP_f", value=f"{delta_P_f_pa/1000.0:,.2f} kPa", delta=f"{delta_P_f_pa/1e5:.3f} bar")

with m_col6:
    st.metric(label="Elevation ΔP_z", value=f"{delta_P_z_pa/1000.0:,.2f} kPa", delta=f"Δz = {delta_z_m:.1f} m")

with m_col7:
    st.metric(label="Total Required ΔP_total", value=f"{delta_P_total_pa/1000.0:,.2f} kPa", delta=f"{delta_P_total_pa/1e5:.3f} bar")

with m_col8:
    st.metric(label="Hydraulic Power", value=f"{p_kw:.2f} kW", delta=f"{p_hp:.2f} HP")


# -----------------------------------------------------------------------------
# DETAILED TABBED RESULTS & VISUALIZATIONS
# -----------------------------------------------------------------------------
tab_table, tab_chart_opt, tab_chart_moody = st.tabs(
    ["📋 Detailed Results Table", "📈 Operating Sensitivity Curve", "🌀 Interactive Moody Diagram"]
)

# --- TAB 1: PANDAS RESULTS TABLE ---
with tab_table:
    st.markdown("### Comprehensive Hydraulic Analysis Data")
    
    table_data = [
        {"Parameter": "Volumetric Flow Rate", "Symbol": "Q", "Value": f"{q_m3s:.6f}", "Unit": "m³/s", "Context": f"Inputs: {input_q_val:.4f} {q_unit}"},
        {"Parameter": "Internal Pipe Diameter", "Symbol": "D", "Value": f"{diameter_m:.4f}", "Unit": "m", "Context": f"Inputs: {input_dia_val:.2f} {dia_unit}"},
        {"Parameter": "Pipe Length", "Symbol": "L", "Value": f"{length_m:.2f}", "Unit": "m", "Context": f"Inputs: {input_len_val:.2f} {len_unit}"},
        {"Parameter": "Pipe Surface Roughness", "Symbol": "ε", "Value": f"{roughness_m*1000.0:.4f}", "Unit": "mm", "Context": f"Relative roughness ε/D = {relative_roughness:.6f}"},
        {"Parameter": "Elevation Difference", "Symbol": "Δz", "Value": f"{delta_z_m:.2f}", "Unit": "m", "Context": "Positive = uphill flow; Negative = downhill flow"},
        {"Parameter": "Fluid Density", "Symbol": "ρ", "Value": f"{rho_kg_m3:.2f}", "Unit": "kg/m³", "Context": f"Selected preset: {fluid_choice}"},
        {"Parameter": "Dynamic Viscosity", "Symbol": "μ", "Value": f"{mu_pas:.6f}", "Unit": "Pa·s", "Context": f"Viscosity = {mu_pas*1000.0:.3f} cP"},
        {"Parameter": "Cross-Sectional Area", "Symbol": "A", "Value": f"{area_m2:.6f}", "Unit": "m²", "Context": "A = πD²/4"},
        {"Parameter": "Average Flow Velocity", "Symbol": "v", "Value": f"{velocity_m_s:.3f}", "Unit": "m/s", "Context": "v = Q/A"},
        {"Parameter": "Reynolds Number", "Symbol": "Re", "Value": f"{reynolds:,.1f}", "Unit": "—", "Context": "Re = ρvD/μ"},
        {"Parameter": "Flow Regime", "Symbol": "—", "Value": flow_regime, "Unit": "—", "Context": "Laminar (<2300), Trans (2300-4000), Turb (≥4000)"},
        {"Parameter": "Darcy Friction Factor", "Symbol": "f", "Value": f"{f_factor:.6f}", "Unit": "—", "Context": method_desc},
        {"Parameter": "Frictional Head Loss", "Symbol": "h_f", "Value": f"{head_loss_m:.3f}", "Unit": "m", "Context": "Darcy-Weisbach: h_f = f(L/D)(v²/2g)"},
        {"Parameter": "Frictional Pressure Loss", "Symbol": "ΔP_f", "Value": f"{delta_P_f_pa/1000.0:.2f}", "Unit": "kPa", "Context": "ΔP_f = ρ g h_f"},
        {"Parameter": "Elevation Pressure Change", "Symbol": "ΔP_z", "Value": f"{delta_P_z_pa/1000.0:.2f}", "Unit": "kPa", "Context": "ΔP_z = ρ g Δz"},
        {"Parameter": "Total Required Pressure Drop", "Symbol": "ΔP_total", "Value": f"{delta_P_total_pa/1000.0:.2f}", "Unit": "kPa", "Context": "ΔP_total = ΔP_f + ΔP_z"},
        {"Parameter": "Hydraulic Power", "Symbol": "P_hyd", "Value": f"{p_kw:.3f}", "Unit": "kW", "Context": f"P_hydraulic = ΔP_total * Q ({p_hp:.2f} HP)"},
    ]

    df_results = pd.DataFrame(table_data)
    st.dataframe(df_results, hide_index=True)

    csv_data = df_results.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 Download Engineering Results Report (CSV)",
        data=csv_data,
        file_name="fluid_flow_calculator_results.csv",
        mime="text/csv",
    )


# --- TAB 2: OPERATING SENSITIVITY CURVE ---
with tab_chart_opt:
    st.markdown("### Pressure Losses vs. Volumetric Flow Rate")
    st.caption("Evaluates system pressure drop response across 10% to 200% of current operating flow rate.")

    # Generate 50 points of Q
    q_range = np.linspace(0.1 * q_m3s, 2.0 * q_m3s, 50)
    dp_f_list = []
    dp_total_list = []
    v_list = []

    for q_i in q_range:
        v_i = q_i / area_m2
        re_i = (rho_kg_m3 * v_i * diameter_m) / mu_pas
        f_i, _, _ = calculate_friction_factor(re_i, relative_roughness)
        hf_i = calculate_head_loss(f_i, length_m, diameter_m, v_i)
        dp_f_i, dp_z_i, dp_tot_i = calculate_pressure_drops(rho_kg_m3, hf_i, delta_z_m)
        dp_f_list.append(dp_f_i / 1000.0)  # kPa
        dp_total_list.append(dp_tot_i / 1000.0)  # kPa
        v_list.append(v_i)

    # Convert x-axis for display in user selected Q unit
    if q_unit == "m³/s":
        q_disp = q_range
    elif q_unit == "L/s":
        q_disp = q_range * 1000.0
    elif q_unit == "m³/h":
        q_disp = q_range * 3600.0
    elif q_unit == "GPM (gal/min)":
        q_disp = q_range / 6.30902e-5
    elif q_unit == "bbl/day":
        q_disp = q_range / 1.84013e-6

    fig_sens = go.Figure()

    fig_sens.add_trace(
        go.Scatter(
            x=q_disp,
            y=dp_f_list,
            mode="lines",
            name="Frictional Loss ΔP_f (kPa)",
            line=dict(color="#2563EB", width=2.5),
        )
    )

    fig_sens.add_trace(
        go.Scatter(
            x=q_disp,
            y=dp_total_list,
            mode="lines",
            name="Total Required ΔP_total (kPa)",
            line=dict(color="#DC2626", width=2.5, dash="dash"),
        )
    )

    # Operating Point Marker
    fig_sens.add_trace(
        go.Scatter(
            x=[input_q_val],
            y=[delta_P_total_pa / 1000.0],
            mode="markers",
            name="Current Operating Point",
            marker=dict(color="#10B981", size=14, symbol="star"),
        )
    )

    fig_sens.update_layout(
        title=f"System Pressure Drop Sensitivity Curve ({fluid_choice})",
        xaxis_title=f"Volumetric Flow Rate Q ({q_unit})",
        yaxis_title="Pressure Drop (kPa)",
        hovermode="x unified",
        template="plotly_white",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01),
        margin=dict(l=40, r=40, t=40, b=40),
    )

    st.plotly_chart(fig_sens)


# --- TAB 3: INTERACTIVE MOODY DIAGRAM ---
with tab_chart_moody:
    st.markdown("### Interactive Moody Diagram Coordinate Mapping")
    st.caption("Standard Moody Chart showing the calculated operating point (Re, f) relative to relative roughness iso-lines.")

    re_arr_lam = np.logspace(2.8, np.log10(2300), 30)
    f_arr_lam = 64.0 / re_arr_lam

    re_arr_turb = np.logspace(np.log10(4000), 8, 100)

    fig_moody = go.Figure()

    # Laminar Line
    fig_moody.add_trace(
        go.Scatter(
            x=re_arr_lam,
            y=f_arr_lam,
            mode="lines",
            name="Laminar (f = 64/Re)",
            line=dict(color="#000000", width=2.5),
        )
    )

    # Iso-roughness curves
    roughness_curves = [0.0, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.03, 0.05]

    for ed in roughness_curves:
        f_curve = []
        for re_val in re_arr_turb:
            f_val, _, _ = calculate_colebrook_white_friction_factor(re_val, ed)
            f_curve.append(f_val)

        label = "Smooth Pipe (e/D=0)" if ed == 0.0 else f"e/D = {ed}"
        fig_moody.add_trace(
            go.Scatter(
                x=re_arr_turb,
                y=f_curve,
                mode="lines",
                name=label,
                line=dict(width=1.2),
            )
        )

    # Current Operating Point Marker
    fig_moody.add_trace(
        go.Scatter(
            x=[reynolds],
            y=[f_factor],
            mode="markers+text",
            name="Current Operating Point",
            text=[f"  Re={reynolds:,.0f}, f={f_factor:.4f}"],
            textposition="top right",
            marker=dict(color="#EF4444", size=14, symbol="diamond"),
        )
    )

    fig_moody.update_layout(
        title="Darcy Friction Factor vs. Reynolds Number (Moody Chart)",
        xaxis_title="Reynolds Number (Re)",
        yaxis_title="Darcy Friction Factor (f)",
        xaxis_type="log",
        yaxis_type="log",
        template="plotly_white",
        margin=dict(l=40, r=40, t=40, b=40),
        legend=dict(font=dict(size=10)),
    )

    st.plotly_chart(fig_moody)

# -----------------------------------------------------------------------------
# FOOTER
# -----------------------------------------------------------------------------
st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: #64748B; font-size: 0.85rem;'>"
    "Fluid Flow Calculator • Petroleum & Mechanical Engineering Hydraulics Tool • Powered by Streamlit"
    "</div>",
    unsafe_allow_html=True,
)
