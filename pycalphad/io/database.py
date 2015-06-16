"""The database module provides support for reading and writing data types
associated with structured thermodynamic/kinetic data.
"""
try:
    set
except NameError:
    from sets import Set as set #pylint: disable=W0622
from tinydb import TinyDB
from tinydb.storages import MemoryStorage

class Database(object): #pylint: disable=R0902
    """
    Structured thermodynamic and/or kinetic data.
    Databases are usually filled by Parsers and read by Models.

    Attributes
    ----------
    elements : list
        List of elements in database.
    species : list
        List of species in database.
    phases : dict
        Phase objects indexed by their system-local name.
    symbols : dict
        SymPy objects indexed by their name (FUNCTIONs in Thermo-Calc).
    references : dict
        Reference objects indexed by their system-local identifier.

    Methods
    -------
    None yet.

    """
    class Phase(object): #pylint: disable=R0903
        """
        Phase in the database.

        Attributes
        ----------
        name : string
            System-local name of the phase.
        constituents : list of lists
            Possible sublattice constituents (elements and/or species).
        sublattices : list
            Site ratios of sublattices.
        model_hints : dict
            Structured "hints" for a Model trying to read this phase.
            Hints for major constituents and typedefs (Thermo-Calc) go here.
        """
        def __init__(self):
            self.name = None
            self.constituents = None
            self.sublattices = []
            self.model_hints = {}
        def __repr__(self):
            return 'Phase({0!r})'.format(self.__dict__)
    def __init__(self, *dbf):
        """
        Construct a Database object.

        Parameters
        ----------
        dbf: file descriptor or raw data, optional
            TDB file to load.

        Examples
        --------
        >>> mydb = Database(open('crfeni_mie.tdb'))
        >>> mydb = Database('crfeni_mie.tdb')
        >>> f = io.StringIO(u'$a complete TDB file as a string\n')
        >>> mydb = Database(f)
        """
        # Should elements be rolled into a special case of species?
        self.elements = set()
        self.species = set()
        self.phases = {}
        self.typedefs = {}
        self._structure_dict = {} # System-local phase names to global IDs
        self._parameters = TinyDB(storage=MemoryStorage)
        self.symbols = {}
        self.references = {}
        # Note: No public typedefs here (from TDB files)
        # Instead we put that information in the model_hint for phases

        if len(dbf) == 1:
            # Let's try to load this database file
            mydb = dbf[0]
            raw_data = ''
            # First, let's try to treat it like it's a file descriptor
            try:
                raw_data = mydb.read()
            except AttributeError:
                # It's not file-like, so it's probably string-like
                # The question is if it's a filename or the whole raw database
                # Solution is to check for newlines
                if mydb.find('\n') == -1:
                    # Single-line; it's probably a filename
                    raw_data = open(mydb).read()
                else:
                    # Newlines found: probably a full database string
                    raw_data = mydb
            # Raw data should be loaded now
            # File type detection (TDB, etc.) would go here
            from pycalphad.io.tdb import tdbread
            tdbread(self, raw_data)
        elif len(dbf) > 1:
            raise ValueError('Invalid number of parameters: '+len(dbf))

    def __str__(self):
        result = 'Elements: {0}\n'.format(sorted(self.elements))
        result += 'Species: {0}\n'.format(sorted(self.species))
        for symbol, info in sorted(self.typedefs.items()):
            result += 'Type Definition \'{0}\': {1}\n'.format(symbol, info)
        for name, phase in sorted(self.phases.items()):
            result += str(phase)+'\n'
        result += '{0} symbols in database\n'.format(len(self.symbols))
        result += '{0} parameters in database\n'.format(len(self._parameters))
        return result

    def add_structure_entry(self, local_name, global_name):
        """
        Define a relation between the system-local name of a phase and a
        "global" identifier. This is used to link crystallographically
        similar phases known by different colloquial names.

        Parameters
        ----------
        local_name : string
            System-local name of the phase.
        global_name : object
            Abstract representation of symbol, e.g., in SymPy format.

        Examples
        --------
        None yet.
        """
        self._structure_dict[local_name] = global_name
    def add_parameter(self, param_type, phase_name, #pylint: disable=R0913
                      constituent_array, param_tag, param_order,
                      param, ref=None):
        """
        Add a parameter.

        Parameters
        ----------
        param_type : str
            Type name of the parameter, e.g., G, L, BMAGN.
        phase_name : string
            Name of the phase.
        constituent_array : list
            Configuration of the sublattices (elements and/or species).
        param_tag : str
            The indetification for mobility parameter, e.g, MQ, MF 
        symbol : object
            Abstract representation of the parameter, e.g., in SymPy format.

        Examples
        --------
        None yet.
        """
        new_parameter = {
            'phase_name': phase_name,
            'constituent_array': constituent_array,
            'parameter_type': param_type,
            'parameter_tag': param_tag,
            'parameter_order': param_order,
            'parameter': param,
            'reference': ref
        }
        param_id = self._parameters.insert(new_parameter)
        return param_id
    def add_phase(self, phase_name, model_hints, sublattices):
        """
        Add a phase.

        Parameters
        ----------
        phase_name : string
            System-local name of the phase.
        model_hints : list
            Structured "hints" for a Model trying to read this phase.
            Hints for major constituents and typedefs (Thermo-Calc) go here.
        sublattices : list
            Site ratios of sublattices.

        Examples
        --------
        None yet.
        """
        new_phase = Database.Phase()
        new_phase.name = phase_name
        new_phase.sublattices = sublattices
        new_phase.model_hints = model_hints
        self.phases[phase_name] = new_phase
    def add_phase_constituents(self, phase_name, constituents):
        """
        Add a phase.

        Parameters
        ----------
        phase_name : string
            System-local name of the phase.
        constituents : list
            Possible phase constituents (elements and/or species).

        Examples
        --------
        None yet.
        """
        try:
            self.phases[phase_name].constituents = constituents
        except KeyError:
            print("Undefined phase "+phase_name)
            raise
    def search(self, query):
        """
        Search for parameters matching the specified query.

        Parameters
        ----------
        query : object
            Structured database query in TinyDB format.

        Examples
        --------
        >>>> from tinydb import where
        >>>> db = Database()
        >>>> eid = db.add_parameter(...) #TODO
        >>>> db.search(where('eid') == eid)
        """
        return self._parameters.search(query)

if __name__ == "__main__":
    pass
