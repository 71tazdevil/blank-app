import streamlit as st
import numpy as np
import astropy as ap
import matplotlib.pyplot as plt
import astroplan

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

ARO = EarthLocation(lat=38.45 * u.deg, lon=-7.54 * u.deg, height=269 * u.m)
time = Time("2024-12-13 20:36:00")


st.title("🎈 My first streamlit app")
st.write(
    "Let's start building! For help and inspiration, head over to [docs.streamlit.io](https://docs.streamlit.io/)."
)

index_messier = 81

cible_coord = SkyCoord.from_name("m" + str(index_messier))



cible_altaz = cible_coord.transform_to(AltAz(obstime=time, location=ARO))
print(f"M33's Altitude = {cible_altaz.alt:.2}")


midi = Time("2024-12-12 17:00:00")
delta_midi = np.linspace(0, 15, 500) * u.hour
frame_date = AltAz(obstime=midi + delta_midi, location=ARO)
cible_altazs_July13night = cible_coord.transform_to(frame_date)

plt.plot(delta_midi, cible_altazs_July13night.alt, color="r", label="M33")
plt.show()








