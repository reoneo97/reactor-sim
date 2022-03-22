import streamlit as st
from simulator.reaction import ptsa_reaction, ReactionModel
import simulator.const as const
from simulator.reactor import IdealCSTR, IdealPFR, IdealPBR
from loguru import logger
import plotly.express as px
import pandas as pd


def conversion_plot(conversion_slc, vol_intervals, time):
    df = pd.DataFrame([vol_intervals, conversion_slc, ]).T
    df.columns = ["Volume (m^3)", "Conversion"]
    fig = px.line(
        df, x="Volume (m^3)", y="Conversion",
        title=f"Conversion Plot for t={time} min")
    final_x, final_y = df.iloc[-1]
    fig.add_annotation(x=final_x, y=final_y,
                       text=f"Output Conversion={final_y:.4f}",
                       showarrow=True,
                       arrowhead=3,
                       xshift=-15)
    st.plotly_chart(fig)


def pbr():
    # Session State Variables
    if 'simulation_done' not in st.session_state:
        st.session_state['simulation_done'] = False
        st.session_state["conversion"] = None
        st.session_state["vol_intervals"] = None

    st.header("Ideal Packed Bed Reactor (PBR) Design")
    st.subheader("[WIP] ⚠️ Not Completed!")
    st.markdown(
        """Ideal PBR design makes several assumptions:

        1. Concentration profile within the PFR is radially uniform
        2. Reactor is perfectly isothermal 
        3. There is no back-mixing within the reactor

        """
    )

    st.subheader("Choosing Design Parameters:")

    temp = st.slider("Temperature (K)", min_value=350.,
                     max_value=480., value=393.15)
    pa_feed_hr = st.slider("Palmitic Acid Feed (kmol//hr)",
                           min_value=1., max_value=100., value=80.31)
    pa_feed = pa_feed_hr/60
    M = st.slider("Molar Ratio (Isopropyl Alcohol/Palmitic Acid)",
                  min_value=1, max_value=5)
    vol = st.slider("Reactor Volume", min_value=0.5,
                    max_value=100., step=0.5, value=35.)

    weight = st.slider("Catalyst Weight", min_value=1000.,
                       max_value=1e10, step=1000.)
    rate_test = st.slider("Rate Coeff", min_value=1, max_value=100, step=1)
    time_interval = 0.05
    time_end = 100.
    st.write(
        f"Using Default Simulation of {time_end} min in {time_interval} min timesteps")
    sim_btn = st.button("Run Simulation")
    if sim_btn:

        model = IdealPBR(pa_feed, M, vol, weight, temp, rate_test=rate_test)
        with st.spinner("Running Simulation"):
            model.run(time_interval, time_end)
        conversion = model.get_conversion()
        st.session_state['simulation_done'] = True
        st.session_state["conversion"] = conversion
        st.session_state["vol_intervals"] = model.get_volume_intervals()
        # st.write(conversion.shape)

    st.markdown("---")
    if st.session_state["simulation_done"]:
        with st.container():
            st.header("Conversion Visualization 📈")
            conversion = st.session_state["conversion"]

            time_slider = st.slider(
                "Time Slider", 0., time_end, value=0., step=time_interval,)

            time_index = min(int(time_slider/time_interval), len(conversion)-1)

            conversion_slc = conversion[time_index, :]
            conversion_plot(
                conversion_slc, st.session_state["vol_intervals"], time_slider)
            # Calculate the Time Interval
            # Plot the Concentration Profile
            # Print out the outlet concentration
            # Find the volume value where

    st.markdown("---")
