import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time, TimeDelta
from astropy import units as u


LCO_LATITUDE = -29.0146 
LCO_LONGITUDE = 289.3074
LCO_ALTITUDE = 2380

def TimeRange(start_time, end_time, time_step):
    """Function to generate an array of time values between start and end with a given time step"""
    start = Time(start_time)
    end = Time(end_time)
    step = TimeDelta(time_step*u.minute)
    return np.arange(start, end, step)
    
    
def plot_star_altitude(
    star_ra: float,
    star_dec: float,
    epoch,
    observer_lat: float,
    observer_lon: float,
    observer_altitude: float,
    start_time: str,
    end_time: str,
    time_step: float,
) -> None:
    """Plot the altitude of a star at a given location and range of dates."""
    
    # Define the observing location
    observing_location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg, height=observer_altitude*u.m)

    # Define the star coordinates with the specified epoch
    star_coords = SkyCoord(ra=star_ra*u.deg, dec=star_dec*u.deg, frame='icrs', equinox=epoch)

    # Generate an array of time values
    times = TimeRange(start_time, end_time, time_step)

    # Convert the time values to the altitude-azimuth frame
    altaz_frame = AltAz(obstime=times, location=observing_location)
    star_altaz = star_coords.transform_to(altaz_frame)

    # Extract the altitude values
    altitudes = star_altaz.alt.degree

    # Convert times to datetime objects
    times_datetime = [t.datetime for t in times]

    # Plot the altitude as a function of time
    plt.plot(times_datetime, altitudes)
    plt.xlabel('Time')
    plt.ylabel('Altitude (degrees)')
    plt.title('Star Altitude Throughout the Night')
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def main():
    """Example usage"""
    star_ra = 197.5  # Right Ascension of the star in degrees
    star_dec = -1.5  # Declination of the star in degrees
    epoch = 'J2000'  # Epoch of the star coordinates
    observer_lat = -29.015972  # Latitude of the observing location in degrees
    observer_lon = -70.69208  # Longitude of the observing location in degrees
    observer_altitude = 2400  # Altitude of the observing location in meters
    start_time = '2023-07-19T22:00:00'  # Start time of the observation (YYYY-MM-DDTHH:MM:SS) UT
    end_time = '2023-07-20T06:00:00'  # End time of the observation (YYYY-MM-DDTHH:MM:SS) UT
    time_step = 10  # Time step between measurements in minutes

    plot_star_altitude(star_ra, star_dec, epoch, observer_lat, observer_lon, observer_altitude, start_time, end_time, time_step)


if __name__ == '__main__':
    main()
