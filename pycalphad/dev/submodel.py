"""
The Model module provides support for using a Database to perform
calculations under specified conditions.
"""
from __future__ import division
import copy
from sympy import log, Add, Mul, Piecewise, Pow, S, Symbol, exp
from tinydb import where
import pycalphad.variables as v
from pycalphad.log import logger
from pycalphad.io.database import Database
from pycalphad.io.utils import print_ccode
try:
    set
except NameError:
    from sets import Set as set #pylint: disable=W0622

class DofError(Exception):
    "Error due to missing degrees of freedom."
    pass

class SubModel(object):
    """
    """
    def __init__(self, dbe, comps, phase, parameters=None):
        possible_comps = set([x.upper() for x in comps])
        print possible_comps
        self.components = set()
        for sublattice in dbe.phases[phase.upper()].constituents:
            self.components |= set(sublattice).intersection(possible_comps)
        logger.debug('Model of %s has components %s', phase, self.components)
        print self.components
        # Verify that this phase is still possible to build
        for sublattice in dbe.phases[phase.upper()].constituents:
            if len(set(sublattice).intersection(self.components)) == 0:
                # None of the components in a sublattice are active
                # We cannot build a model of this phase
                raise DofError(
                    '{0}: Sublattice {1} of {2} has no components in {3}' \
                    .format(phase.upper(), sublattice,
                            dbe.phases[phase.upper()].constituents,
                            self.components))
        symbols = dict([(Symbol(s), val) for s, val in dbe.symbols.items()])
        if parameters is not None:
            symbols.update([(Symbol(s), val) for s, val in parameters.items()])
        # Need to do more substitutions to catch symbols that are functions
        # of other symbols
        for name, value in symbols.items():
            try:
                symbols[name] = value.xreplace(symbols)
            except AttributeError:
                # Can't use xreplace on a float
                pass
        for name, value in symbols.items():
            try:
                symbols[name] = value.xreplace(symbols)
            except AttributeError:
                # Can't use xreplace on a float
                pass
        self.ast = self.build_phase(dbe, phase.upper(), symbols, dbe.search)
    def _purity_test(self, constituent_array):
        """
        Check if constituent array only has one species in its array
        This species must also be an active species
        """
        for sublattice in constituent_array:
            if len(sublattice) != 1:
                return False
            if (sublattice[0] not in self.components) and \
                (sublattice[0] != '*'):
                return False
        return True
    def _array_validity(self, constituent_array):
        """
        Check that the current array contains only active species.
        """
        for sublattice in constituent_array:
            valid = set(sublattice).issubset(self.components) \
                or sublattice[0] == '*'
            if not valid:
                return False
        return True
    def _interaction_test(self, constituent_array):
        """
        Check if constituent array has more than one active species in
        its array for at least one sublattice.
        """
        result = False
        for sublattice in constituent_array:
            # check if all elements involved are also active
            valid = set(sublattice).issubset(self.components) \
                or sublattice[0] == '*'
            if len(sublattice) > 1 and valid:
                result = True
            if not valid:
                result = False
                break
        return result

    def _site_ratio_normalization(self, phase):
        """
        Calculates the normalization factor based on the number of sites
        in each sublattice.
        """
        site_ratio_normalization = S.Zero
        # Normalize by the sum of site ratios times a factor
        # related to the site fraction of vacancies
        for idx, sublattice in enumerate(phase.constituents):
            if ('VA' in set(sublattice)) and ('VA' in self.components):
                site_ratio_normalization += phase.sublattices[idx] * \
                    (1 - v.SiteFraction(phase.name, idx, 'VA'))
            else:
                site_ratio_normalization += phase.sublattices[idx]
        return site_ratio_normalization

    @staticmethod
    def _Muggianu_correction_dict(comps): #pylint: disable=C0103
        """
        Replace y_i -> y_i + (1 - sum(y involved in parameter)) / m,
        where m is the arity of the interaction parameter.
        Returns a dict converting the list of Symbols (comps) to this.
        m is assumed equal to the length of comps.

        When incorporating binary, ternary or n-ary interaction parameters
        into systems with more than n components, the sum of site fractions
        involved in the interaction parameter may no longer be unity. This
        breaks the symmetry of the parameter. The solution suggested by
        Muggianu, 1975, is to renormalize the site fractions by replacing them
        with a term that will sum to unity even in higher-order systems.
        There are other solutions that involve retaining the asymmetry for
        physical reasons, but this solution works well for components that
        are physically similar.

        This procedure is based on an analysis by Hillert, 1980,
        published in the Calphad journal.
        """
        arity = len(comps)
        return_dict = {}
        correction_term = (S.One - Add(*comps)) / arity
        for comp in comps:
            return_dict[comp] = comp + correction_term
        return return_dict
    def _redlich_kister_sum(self, phase, symbols, param_type, param_search):
        """
        Construct parameter in Redlich-Kister polynomial basis, using
        the Muggianu ternary parameter extension.
        """
        rk_terms = []
        param_query = (
            (where('phase_name') == phase.name) & \
            (where('parameter_type') == param_type) & \
            (where('constituent_array').test(self._array_validity))
        )
        # search for desired parameters
        params = param_search(param_query)
        for param in params:
            # iterate over every sublattice
            mixing_term = S.One
            for subl_index, comps in enumerate(param['constituent_array']):
                comp_symbols = None
                # convert strings to symbols
                if comps[0] == '*':
                    # Handle wildcards in constituent array
                    comp_symbols = \
                        [
                            v.SiteFraction(phase.name, subl_index, comp)
                            for comp in set(phase.constituents[subl_index])\
                                .intersection(self.components)
                        ]
                    mixing_term *= Add(*comp_symbols)
                else:
                    comp_symbols = \
                        [
                            v.SiteFraction(phase.name, subl_index, comp)
                            for comp in comps
                        ]
                    mixing_term *= Mul(*comp_symbols)
                # is this a higher-order interaction parameter?
                if len(comps) == 2 and param['parameter_order'] > 0:
                    # interacting sublattice, add the interaction polynomial
                    redlich_kister_poly = Pow(comp_symbols[0] - \
                        comp_symbols[1], param['parameter_order'])
                    mixing_term *= redlich_kister_poly
                if len(comps) == 3:
                    # 'parameter_order' is an index to a variable when
                    # we are in the ternary interaction parameter case

                    # NOTE: The commercial software packages seem to have
                    # a "feature" where, if only the zeroth
                    # parameter_order term of a ternary parameter is specified,
                    # the other two terms are automatically generated in order
                    # to make the parameter symmetric.
                    # In other words, specifying only this parameter:
                    # PARAMETER G(FCC_A1,AL,CR,NI;0) 298.15  +30300; 6000 N !
                    # Actually implies:
                    # PARAMETER G(FCC_A1,AL,CR,NI;0) 298.15  +30300; 6000 N !
                    # PARAMETER G(FCC_A1,AL,CR,NI;1) 298.15  +30300; 6000 N !
                    # PARAMETER G(FCC_A1,AL,CR,NI;2) 298.15  +30300; 6000 N !
                    #
                    # If either 1 or 2 is specified, no implicit parameters are
                    # generated.
                    # We need to handle this case.
                    if param['parameter_order'] == 0:
                        # are _any_ of the other parameter_orders specified?
                        ternary_param_query = (
                            (where('phase_name') == param['phase_name']) & \
                            (where('parameter_type') == \
                                param['parameter_type']) & \
                            (where('constituent_array') == \
                                param['constituent_array'])
                        )
                        other_tern_params = param_search(ternary_param_query)
                        if len(other_tern_params) == 1 and \
                            other_tern_params[0] == param:
                            # only the current parameter is specified
                            # We need to generate the other two parameters.
                            order_one = copy.deepcopy(param)
                            order_one['parameter_order'] = 1
                            order_two = copy.deepcopy(param)
                            order_two['parameter_order'] = 2
                            # Add these parameters to our iteration.
                            params.extend((order_one, order_two))
                    # Include variable indicated by parameter order index
                    # Perform Muggianu adjustment to site fractions
                    mixing_term *= comp_symbols[param['parameter_order']].subs(
                        self._Muggianu_correction_dict(comp_symbols),
                        simultaneous=True)
            rk_terms.append(mixing_term * param['parameter'].xreplace(symbols))
        return Add(*rk_terms)
    def build_phase(self, dbe, phase_name, symbols, param_search):
        """
        Apply phase's model hints to build a master SymPy object.
        """
        phase = dbe.phases[phase_name]
        mobility = {}
        # First, build the reference energy term
        for eachcom in self.components:

            mf = self.pure_contribution(
                phase, symbols, eachcom, "MF", param_search
            )
            mf += self.excess_mixing_energy(
                phase, symbols, eachcom, "MF", param_search
            )
            mq = self.pure_contribution(
                phase, symbols, eachcom, "MQ", param_search
            )
            mq += self.excess_mixing_energy(
                phase, symbols, eachcom, "MQ", param_search
            )
            mf = mf == S.Zero and S.One or mf
            mobility[eachcom] = exp(mf/v.R/v.T)*exp(mq/v.R/v.T)/v.R/v.T
        return mobility
    def pure_contribution(self, phase, symbols, param_tag,
        param_type, param_search):
        """
        Returns the weighted average of the endmember energies
        in symbolic form.
        """
        pure_energy_term = S.Zero
        site_ratio_normalization = self._site_ratio_normalization(phase)

        pure_param_query = (
            (where('phase_name') == phase.name) & \
            (where('parameter_order') == 0) & \
            (where('parameter_type') == param_type) & \
            (where('parameter_tag') == param_tag) & \
            (where('constituent_array').test(self._purity_test))
        )

        pure_params = param_search(pure_param_query)

        for param in pure_params:
            site_fraction_product = S.One
            for subl_index, comp in enumerate(param['constituent_array']):
                # We know that comp has one entry, by filtering
                if comp[0] == '*':
                    # Handle wildcards in constituent array
                    comp_symbols = \
                        [
                            v.SiteFraction(phase.name, subl_index, compx)
                            for compx in set(phase.constituents[subl_index])\
                                .intersection(self.components)
                        ]
                    site_fraction_product *= Add(*comp_symbols)
                else:
                    if comp[0] not in self.components:
                        # not a valid parameter for these components
                        # zero out this contribution
                        site_fraction_product = S.Zero
                    else:
                        comp_symbol = \
                            v.SiteFraction(phase.name, subl_index, comp[0])
                        site_fraction_product *= comp_symbol
            pure_energy_term += (
                site_fraction_product * param['parameter'].xreplace(symbols)
            ) / site_ratio_normalization
        return pure_energy_term
    def excess_mixing_energy(self, phase, symbols, param_tag,
        param_type, param_search):
        """
        Build the binary, ternary and higher order interaction term
        Here we use Redlich-Kister polynomial basis by default
        Here we use the Muggianu ternary extension by default
        Replace y_i -> y_i + (1 - sum(y involved in parameter)) / m,
        where m is the arity of the interaction parameter
        """
        excess_mixing_terms = []
        # Normalize site ratios
        site_ratio_normalization = self._site_ratio_normalization(phase)
        site_ratios = phase.sublattices
        site_ratios = [c/site_ratio_normalization for c in site_ratios]

        interaction_param_query = (
            (where('phase_name') == phase.name) & \
            (where('parameter_type') == param_type) & \
            (where('param_tag') == param_tag) & \
            (where('constituent_array').test(self._interaction_test))
        )
        # search for desired parameters
        interaction_params = param_search(interaction_param_query)

        for param in interaction_params:
            # iterate over every sublattice
            mixing_term = S.One
            for subl_index, comps in enumerate(param['constituent_array']):
                comp_symbols = None
                # convert strings to symbols
                if comps[0] == '*':
                    # Handle wildcards in constituent array
                    comp_symbols = \
                        [
                            v.SiteFraction(phase.name, subl_index, comp)
                            for comp in set(phase.constituents[subl_index])\
                                .intersection(self.components)
                        ]
                    mixing_term *= Add(*comp_symbols)
                else:
                    # Order of elements in comps matters here!
                    # This means we can't call set(comps)
                    # No need to check set intersection here anyway because
                    # only valid parameters will be returned by our query
                    comp_symbols = \
                        [
                            v.SiteFraction(phase.name, subl_index, comp)
                            for comp in comps
                        ]
                    mixing_term *= Mul(*comp_symbols)
                # is this a higher-order interaction parameter?
                if len(comps) == 2 and param['parameter_order'] > 0:
                    # interacting sublattice, add the interaction polynomial
                    redlich_kister_poly = Pow(comp_symbols[0] - \
                        comp_symbols[1], param['parameter_order'])
                    mixing_term *= redlich_kister_poly
                if len(comps) == 3:
                    # 'parameter_order' is an index to a variable when
                    # we are in the ternary interaction parameter case

                    # NOTE: The commercial software packages seem to have
                    # a "feature" where, if only the zeroth
                    # parameter_order term of a ternary parameter is specified,
                    # the other two terms are automatically generated in order
                    # to make the parameter symmetric.
                    # In other words, specifying only this parameter:
                    # PARAMETER G(FCC_A1,AL,CR,NI;0) 298.15  +30300; 6000 N !
                    # Actually implies:
                    # PARAMETER G(FCC_A1,AL,CR,NI;0) 298.15  +30300; 6000 N !
                    # PARAMETER G(FCC_A1,AL,CR,NI;1) 298.15  +30300; 6000 N !
                    # PARAMETER G(FCC_A1,AL,CR,NI;2) 298.15  +30300; 6000 N !
                    #
                    # If either 1 or 2 is specified, no implicit parameters are
                    # generated.
                    # We need to handle this case.
                    if param['parameter_order'] == 0:
                        # are _any_ of the other parameter_orders specified?
                        ternary_param_query = (
                            (where('phase_name') == param['phase_name']) & \
                            (where('parameter_type') == \
                                param['parameter_type']) & \
                            (where('constituent_array') == \
                                param['constituent_array'])
                        )
                        other_tern_params = param_search(ternary_param_query)
                        if len(other_tern_params) == 1 and \
                            other_tern_params[0] == param:
                            # only the current parameter is specified
                            # We need to generate the other two parameters.
                            order_one = copy.deepcopy(param)
                            order_one['parameter_order'] = 1
                            order_two = copy.deepcopy(param)
                            order_two['parameter_order'] = 2
                            # Add these parameters to our iteration.
                            interaction_params.extend((order_one, order_two))
                    # Include variable indicated by parameter order index
                    # Perform Muggianu adjustment to site fractions
                    mixing_term *= comp_symbols[param['parameter_order']].subs(
                        self._Muggianu_correction_dict(comp_symbols),
                        simultaneous=True)
                if len(comps) > 3:
                    raise ValueError('Higher-order interactions (n>3) are \
                        not yet supported')
            excess_mixing_terms.append(mixing_term * \
                param['parameter'].xreplace(symbols))
        return Add(*excess_mixing_terms)/site_ratio_normalization
if __name__ == "__main__":
    import submodel_plugin as plugin

    DBF = Database(plugin.TDB_TEST_STRING)
    submodel = SubModel(DBF, ["AL", "CU", "SI", "VA"], "FCC_A1")
    print submodel.ast
