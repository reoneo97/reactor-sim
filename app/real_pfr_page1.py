import matplotlib.pyplot as plt
import streamlit as st
from simulator.reaction import ptsa_reaction, ReactionModel
import simulator.const as const
from simulator.reactor import RealPFR
from simulator.metrics import total_cost, get_tac
from loguru import logger
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import math
import numpy as np


def conversion_mean_plot(conversion_slc, vol_weights, z_axis, r_axis, time):
    normalised = conversion_slc@vol_weights
    df = pd.DataFrame([z_axis, normalised]).T
    # norm_temp = temp_slc*vol_weights
    df.columns = ["Length/Z-Axis (m)", "Conversion"]
    fig = px.line(
        df, x=df.columns[0], y=df.columns[1],
        title=f"Mean Conversion Plot for t={time} min")
    final_x, final_y = df.iloc[-1]
    fig.add_annotation(x=final_x, y=final_y,
                       text=f"Output Conversion={final_y:.4f}",
                       showarrow=True,
                       arrowhead=3,
                       xshift=-15)
    st.plotly_chart(fig)


def temperature_mean_plot(temp_slc, vol_weights, z_axis, r_axis, time):
    norm_T = temp_slc@vol_weights
    df = pd.DataFrame([z_axis, norm_T]).T
    # norm_temp = temp_slc*vol_weights
    df.columns = ["Length/Z-Axis (m)", "Temperature"]
    fig = px.line(
        df, x=df.columns[0], y=df.columns[1],
        title=f"Mean Temperature Plot for t={time} min")
    final_x, final_y = df.iloc[-1]
    fig.add_annotation(x=final_x, y=final_y,
                       text=f"Output Temperature={final_y:.4f}",
                       showarrow=True,
                       arrowhead=3,
                       xshift=-15)
    st.plotly_chart(fig)


def conversion_contour(conversion_slc, z_axis, r_axis):
    # Plotting is X x Y x Z which will translate to L x R x T
    contour_obj = go.Contour(
        z=conversion_slc.T, x=z_axis, y=r_axis
    )
    fig = go.Figure(contour_obj)
    fig.layout.xaxis.title = "Length/Z-axis Profile (m)"
    fig.layout.yaxis.title = "Radius Profile (m)"
    fig.layout.title = "<b>Conversion Contour</b>"
    st.plotly_chart(fig)


def temperature_contour(temp_slc, z_axis, r_axis):
    # Plotting is X x Y x Z which will translate to L x R x T
    contour_obj = go.Contour(
        z=temp_slc.T, x=z_axis, y=r_axis
    )
    fig = go.Figure(contour_obj)
    fig.layout.xaxis.title = "Length/Z-axis Profile (m)"
    fig.layout.yaxis.title = "Radius Profile (m)"
    fig.layout.title = "<b>Temperature Contour</b>"
    st.plotly_chart(fig)


def heater_requirement_plot(time, heat_flow):
    df = pd.DataFrame([time, heat_flow, ]).T
    df.columns = ["Time (min)", "Heating Requirement (kJ/min)"]
    fig = px.line(
        df, x="Time (min)", y="Heating Requirement (kJ/min)",
        title=f"Heating Requirement Plot")
    st.plotly_chart(fig)


def real_pfr_iso():

    # Session State Variables
    if 'simulation_done' not in st.session_state:
        st.session_state['simulation_done'] = False
        st.session_state["simulation_data"] = None
        st.session_state["vol_intervals"] = None
        st.session_state["sim_time_end"] = None
        st.session_state["heater_flow"] = None
        st.session_state["model"] = None
        st.session_state["costing"] = None
        st.session_state["TAC"] = None

    st.header("Real PFR Design")
    # st.subheader("[WIP] ⚠️ Not Completed!")
    st.markdown(
        """Real PFR design makes several assumptions:

        1. Velocity of Fluid is Radially uniform
        2. There is no back-mixing within the reactor
        3. Heat transfer in radial direction is only limited to conduction. 
        4. Energy loss from heating fluid to the surroundings is negligible

        """
    )

    st.subheader("Choosing Design Parameters:")

    col1, _, col2 = st.columns([8, 1, 8])
    with col1:
        st.markdown("#### 🌡 Reactor Conditions:")
        pa_feed_hr = st.slider("Palmitic Acid Feed (kmol//hr)",
                               min_value=1., max_value=100., value=80.019)
        M = st.slider("Molar Ratio (Isopropyl Alcohol/Palmitic Acid)",
                      min_value=1., max_value=5., value=5., step=0.01)
        feed_temp = st.slider("Feed Temperature (K)",
                              min_value=293., max_value=428.15, value=393., step=0.01)
        heat_temp = st.slider("Heater Temperature (K)",
                              min_value=373., max_value=523., value=433.15)
        n_reactor = st.slider("Number of Parallel Reactors",
                              min_value=1, max_value=8, step=1)
        pa_feed = (pa_feed_hr/60)/n_reactor
        # splitting the pa feed into n parallel reactors
    with col2:
        st.markdown("#### ⚗️ Reactor Dimensions")
        L = st.slider("Reactor Length", min_value=0.1,
                      max_value=25., step=0.1, value=10.)
        DIAMETER = st.slider("Reactor Diameter", min_value=0.4,
                             max_value=6., step=0.01, value=2.)
        R = DIAMETER/2
        reactor_vol = L*math.pi*R*R
        reactor_area = L*math.pi*DIAMETER
        st.write(
            f"Volume: {reactor_vol:.3f} m³ | Surface Area: {reactor_vol:.3f} m²")
        st.markdown("#### Simulation Parameters")
        space_interval = 26
        time_interval = st.slider(
            "Time Step for Simulation", min_value=0.2, max_value=0.5, step=0.05, value=0.2)
        st.info("If the simulation is having issues, reduce the timestep. However this will increase simulation time")
        time_end = st.slider(
            "Simulation Time End", min_value=100., max_value=500., step=50.
        )

    st.write(
        f"Performing Simulation with Time Interval of {time_interval} min")
    sim_btn = st.button("Run Simulation")
    if sim_btn:
        model = RealPFR(pa_feed, M, L, R, feed_temp, heat_temp, space_interval)
        with st.spinner("Running Simulation - Simulation for Real PFR can take up to 5 minutes."):
            sim_data = model.run(time_interval, time_end)

            st.session_state["simulation_done"] = True
            st.session_state["simulation_data"] = sim_data
            st.session_state["sim_time_end"] = time_end
            ls, rs = model.get_dimensions()
            vol_weights = model.get_radial_volumes()
            st.session_state["z_axis"] = ls
            st.session_state["r_axis"] = rs
            st.session_state["vol_weights"] = vol_weights
            st.session_state["heater_flow"] = model.get_heater_total() / \
                time_interval  # Normalise to get in terms of kJ/min
            st.session_state["model"] = model
            cost_dict = total_cost(model, time_interval, n_reactor)

            st.session_state["costing"] = cost_dict
            # st.session_state["TAC"] = get_tac(model, time_interval, n_reactor)

    st.markdown("---")
    cost_container = st.container()
    plot_container = st.container()

    with plot_container:
        if st.session_state["simulation_done"]:
            model = st.session_state["model"]
            st.header("Reactor Output Information")
            simulation_data = st.session_state["simulation_data"]
            vol_weights = st.session_state["vol_weights"]
            output_slc = simulation_data[len(simulation_data)-1, -1, :, :]
            output_temp_slc = output_slc[:, 0]
            output_conc_slc = output_slc[:, 1:6]

            output_temp = output_temp_slc@vol_weights
            output_concs = output_conc_slc.T@vol_weights
            output_conversion = output_slc[:, -1]@vol_weights
            output_df = pd.DataFrame(
                [output_concs],
                columns=[
                    "Palmitic Acid Concentration", "IPA Concentration",
                    "IPP Concentration", "Water Concentration",
                    "Linoleic Acid Concentration"])
            st.write(output_df)
            info = f"""
            Output Temperature: {output_temp:.2f}K \n
            Output Conversion: {output_conversion:.3f} \n
            Fluid Velocity: {model.velocity:.3f} m/min \n
            Pressure Drop: {model.pressure_drop:.3f} kPa \n
            Volumetric Flow Rate: {model.flow_rate:.3f} m3/min \n
            Jacket Heat Flow Rate: {st.session_state['heater_flow'][-1]*60:.3f} kJ/hr \n
            """
            st.write(info)

            st.header("Conversion Visualization 📈")
            st.markdown("---")
            ts = np.arange(0, st.session_state["sim_time_end"], time_interval)
            time_slider = st.slider(
                "Time Slider", 0., st.session_state["sim_time_end"], value=0., step=time_interval,)

            time_index = min(int(time_slider/time_interval),
                             len(simulation_data)-1)

            data_slc = simulation_data[time_index, :]
            plot_col1, _, plot_col2 = st.columns([8, 1, 8])
            z_axis = st.session_state["z_axis"]
            r_axis = st.session_state["r_axis"]
            vol_weights = st.session_state["vol_weights"]
            with plot_col1:
                st.markdown("### **<ins>Temperature Plots</ins>**",
                            unsafe_allow_html=True)
                temp_slc = data_slc[:, :, 0]
                st.markdown("#### Temperature Contour")
                temperature_contour(
                    temp_slc, z_axis, r_axis
                )
                st.markdown("#### Average Temperature Line Plot")
                temperature_mean_plot(
                    temp_slc, vol_weights, z_axis, r_axis, time_slider)
                heater_requirement_plot(ts, st.session_state["heater_flow"])
            with plot_col2:
                st.markdown("### **<ins>Conversion Plots</ins>**",
                            unsafe_allow_html=True)
                conv_slc = data_slc[:, :, -1]
                st.markdown("#### Conversion Contour")
                conversion_contour(
                    conv_slc,
                    st.session_state["z_axis"], st.session_state["r_axis"]
                )
                st.markdown("#### Average Conversion Line Plot")
                conversion_mean_plot(
                    conv_slc, vol_weights, z_axis, r_axis, time_slider)

        if st.session_state["simulation_done"]:
            st.write(pd.DataFrame([st.session_state["costing"]]).T)
