from .pbr_page import pbr
from .pfr_page import pfr
from .real_pfr_page1 import real_pfr_iso
from .optim_pfr import optim

pages = {
    "Ideal PFR (p-TSA)": pfr,
    "Ideal PBR (ZnA/SG)": pbr,
    "Real PFR - Isothermal": real_pfr_iso,
    "Optimized PFR": optim
}
