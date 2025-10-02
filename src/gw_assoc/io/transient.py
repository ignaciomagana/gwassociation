# src/gw_assoc/io/transient.py
class Transient:
    def __init__(self, ra=None, dec=None, z=None, time=None):
        self.ra = ra
        self.dec = dec
        self.z = z
        self.time = time

    def __repr__(self):
        return f"<Transient ra={self.ra}, dec={self.dec}, z={self.z}, time={self.time}>"

