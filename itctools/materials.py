from .itctools import ureg


class Solvent(object):

    """
    A Solvent object represents a liquid that may be pipetted, and in which compounds may be dissolved.

    """

    def __init__(self, name, density=None):
        """
        Parameters
        ----------
        name : str
           The name of the solvent to create.
        density : pint Quantity with units compatible with grams/milliliter, optional, default=None
           The density of the solvent.

        Examples
        --------

        Register a solvent.

        >>> water = Solvent('water', density=0.9970479*ureg.gram/ureg.centimeter**3)

        Register a solvent with density information.

        >>> dmso = Solvent('dmso', density=1.1004*ureg.grams/ureg.centimeter**3)

        """
        self.name = name
        self.density = density


class Compound(object):

    """
    A Compound object represents a compound that can be dissolved in a solvent.

    """

    def __init__(self, name, molecular_weight=None, purity=1.0):
        """
        Parameters
        ----------
        name : str
           The name of the compound to create.
        molecular_weight : pint Quantity with units compatible with grams/mole, optional, default=None
           The molecular weight of the compound.
        purity : float, optional, default=1.0
           The mass purity used for computing actual quantity of compound.

        Examples
        --------

        Register a compound.

        >>> nacl = Compound('sodium chloride')

        Register a compound with molecular weight.

        >>> imatinib = Compound('imatinib mesylate', molecular_weight=589.7*ureg.grams/ureg.mole)

        Use a non-unit purity.

        >>> compound1 = Compound('compound1', molecular_weight=209.12*ureg.grams/ureg.mole, purity=0.975)

        """
        self.name = name
        self.molecular_weight = molecular_weight
        self.purity = purity


class PureLiquid(Compound):

    """A PureLiquid describes a pure liquid that can be part of a mixture of liquids."""

    def __init__(self, name, density, molecular_weight, purity=1.0):
        """
        name : str
            name of the liquid
        density : pint Quantity with units compatible with grams/milliliter
            density of the pure liquid
        molecular weight : pint Quantity with units compatible with grams/mole
            molecular weight of pure liquid
        purity : float, optional, default = 1.0
            fraction of liquid that is pure
        """
        super(
            PureLiquid,
            self).__init__(
            name,
            molecular_weight=molecular_weight,
            purity=purity)
        self.density = density


class SimpleSolution(Solvent):

    """
    A SimpleSolution object represents a solution containing one compound and one solvent.

    The solution is assumed to be ideal, with the same volume as that of the solvent.

    """

    def __init__(
            self,
            compound,
            compound_mass,
            solvent,
            solvent_mass,
            location):
        """
        compound : Compound
           The compound added to the solution.
        compound_mass : pint Quantity compatible with grams
           The mass of compound added to the solution.
        solvent : Solvent
           The solvent used for the solution.
        solvent_mass : pint Quantity compatible with grams
           The mass of solvent used for the solution.
        location : PipettingLocation
           The pipetting location holding the solution.

        Examples
        --------

        Create a simple salt solution.

        >>> salt = Compound('sodium chloride', molecular_weight=58.44277*ureg.grams/ureg.mole)
        >>> water = Solvent('water', density=0.9970479*ureg.grams/ureg.centimeter**3)
        >>> location = PipettingLocation('BufferTrough', 'Trough 100ml', 1)
        >>> solution = SimpleSolution(compound=salt, compound_mass=1.0*ureg.milligrams, solvent=water, solvent_mass=10.0*ureg.grams, location=location)

        TODO
        ----
        * Allow specification of quantity of compound and solvent in various ways (mass, moles, volume) with automated conversions.

        """

        self.compound = compound
        self.compound_mass = compound_mass
        self.solvent = solvent
        self.solvent_mass = solvent_mass

        self.name = compound.name

        # Compute total solution mass.
        self.solution_mass = self.compound_mass + self.solvent_mass

        # Assume solution is ideal; that density and volume is same as solvent.
        self.density = solvent.density
        self.volume = solvent_mass / solvent.density

        # Compute number of moles of compound.
        # number of moles of compound
        self.compound_moles = compound_mass / compound.molecular_weight * compound.purity

        # Compute molarity.
        self.concentration = self.compound_moles / self.volume

        # Store location.
        self.location = location


class SimpleMixture(Solvent):

    """
    A SimpleMixture object represents a solution containing a mixture of various solvents.

    The solution is assumed to be ideal, with the same volume as that of the solvent.

    """

    def __init__(
            self,
            components=list(),
            molefractions=list(),
            locations=list(),
            normalize_fractions=False):
        """
        components : list of PureLiquid
            components of the mixture
        molefractions : list of float
            mole fraction per component
        locations : list of PipettingLocation
            The pipetting location holding the pure liquids
        normalize_fractions : bool, optional, default = False
            Normalize any mole fractions to form a total of 1.
        """
        self.components = components
        self.molefractions = molefractions
        self.locations = locations
        # Consistency checks

        # Input length
        if not len(components) == len(molefractions) == len(locations):
            raise ValueError("Input lists do not have same length!")

        # Ensure total mole fraction equals 1
        if normalize_fractions:
            total = sum(self.molefractions)
            self.molefractions = map(lambda x: x / total, self.molefractions)
        else:
            # Check if mole fraction is 1 within arbitrary precision
            if abs(1.0 - sum(self.molefractions)) > 0.0001:
                raise ValueError("Total mole fractions out of bounds!")

        # Molar volumes are calculated so we can get the volume fractions
        self.molarvolumes = [component.molecular_weight / component.density for component in self.components]

        self.massfractions = self._calculate_mass_fractions()

        self.volumefractions = self._calculate_volume_fractions()


    def _calculate_mass_fractions(self):
        # Average molecular weight  = sum(mole fraction * mole weight) over all components of the mixture
        normalizing_mass = 0 * ureg.grams / ureg.mole
        for i, component in enumerate(self.components):
            normalizing_mass += component.molecular_weight * self.molefractions[i]

        # Mass fraction per component = mole fraction * molecular weight / average molecular weight
        massfractions = [
            self.molefractions[i] *
            component.molecular_weight /
            normalizing_mass for i, component in enumerate(self.components)]

        return massfractions

    def _calculate_volume_fractions(self):
        # Average Molar volume = sum( mole fraction * molar volume) over all components of the mixture
        normalizing_volume = 0 * ureg.liter / ureg.mole
        for i, volume in enumerate(self.molarvolumes):
            normalizing_volume += self.molefractions[i] * volume

        # Volume fraction per component = mole fraction * molar volume / average molar volume
        volumefractions = [
            self.molefractions[i] *
            volume /
            normalizing_volume for i, volume in enumerate(self.molarvolumes)]
        return volumefractions


    def __str__(self):
        """Represent a mixture by its composition."""
        return "<%s: %s>" % (self.__class__, self.describe())

    def describe(self):
        """Give a description of the mixture composition."""
        composition = str()
        for n, comp in enumerate(self.components):
            if self.molefractions[n] > 0.0:
                composition += comp.name
                composition += " %.2f" % self.molefractions[n]
                composition += "; "
        return composition


# MAIN AND TESTS
if __name__ == '__main__':
    import doctest
    doctest.testmod()
