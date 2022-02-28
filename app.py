import streamlit as st
# from simulator.reaction import ptsa_reaction, ReactionModel
import simulator.const as const
from simulator.reactor import IdealCSTR, IdealPFR
from loguru import logger
import plotly.express as px
import pandas as pd
from app import pages

st.set_page_config(layout="wide")
with st.sidebar:

    st.header("Page Navigation")
    pg = st.radio("Page Selection",
                  options=[
                      "Real PFR",
                      "Ideal PFR (p-TSA)",
                      "Ideal PBR (ZnA/SG)",
                  ])


def project_header():
    st.title("Team 25 FYDP - Section 7 Simulator 🧪")
    st.markdown("*Neo Zheng Jie Reo (A0183286W)*")


project_header()
st.markdown("---")
pages[pg]()
