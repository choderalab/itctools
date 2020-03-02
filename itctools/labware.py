
class Labware(object):
    """Define basic Tecan labware."""

    def __init__(self, RackLabel, RackType, RackID=None):
        """
        Parameters
        ----------
        RackLabel : str
           Tecan labware name on deck (e.g. 'SourcePlate', 'Buffer', 'Water')
        RackType : str
           Tecan labeware type (e.g. 'Trough 100ml', 'ITC Plate')
        RackID : str, optional
           Tecan barcode.

        """
        self.RackLabel = RackLabel
        self.RackType = RackType
        self.RackID = RackID


class PipettingLocation(object):
    """Defines a position the Tecan can pipette."""

    def __init__(self, RackLabel, RackType, Position):
        """
        Parameters
        ----------
        RackLabel : str
           Tecan labware name on deck (e.g. 'SourcePlate', 'Buffer', 'Water')
        RackType : str
           Tecan labeware type (e.g. 'Trough 100ml', 'ITC Plate')
        Position :
        TODO finish documentation
        """
        # Information for Tecan LiHa.
        self.RackLabel = RackLabel
        self.RackType = RackType
        self.Position = Position
