/* 
 * Copyright (C) 2015 Gereon Kaiping <gereon.kaiping@soton.ac.uk>
 *
 * This file is part of the BEAST2 package parameterclone.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package parameterclone.splitandmerge;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.Binomial;
import beast.util.Randomizer;

@Description("Randomly split a group of parameters in two")
@Citation("Pagel, M., Meade, A., 2006."
		+" Bayesian Analysis of Correlated Evolution of Discrete Characters"
		+" by Reversible-Jump Markov Chain Monte Carlo."
		+" The American Naturalist 167, 808--825. doi:10.1086/503444")
public class SplitOperator extends Operator {
	// Inputs that are changed by the operator
	public Input<List<RealParameter>> parametersInput = new Input<List<RealParameter>>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			new ArrayList<>(), Validate.REQUIRED);
	public Input<List<IntegerParameter>> groupingsInput = new Input<List<IntegerParameter>>(
			"groupings", "parameter selection indices",
			new ArrayList<>(), Validate.REQUIRED);

	Integer k;
	Integer maxIndex;

	@Override
	public void initAndValidate() throws Exception {
		k = groupingsInput.get().size();
		maxIndex = parametersInput.get().size();
		for (IntegerParameter group : groupingsInput.get()) {
			if (group.getValue() >= maxIndex) {
				throw new Exception(
						"All entries in groupings must be valid indices of parameters");
			}
		}
	}

    /**
     * Change the parameter and return the log of the Hastings ratio.
     * Split a class of joined parameters in two.
     */
	@Override
	public double proposal() {
		// Find the composition of groups, in particular which ones can be split.
		HashMap<Integer, ArrayList<Integer>> groups = new HashMap<Integer, ArrayList<Integer>>();
		ArrayList<Integer> groupsOfSizeAtLeastTwo = new ArrayList<Integer>();
		for (Integer index = 0; index <= groupingsInput.get().size(); ++index) {
			Integer value = groupingsInput.get().get(index).getValue();
			if (groups.get(value) != null) {
				if (!groupsOfSizeAtLeastTwo.contains(value)) {
					groupsOfSizeAtLeastTwo.add(value);
				}
			} else {
				ArrayList<Integer> newGroup = new ArrayList<Integer>();
				groups.put(value, newGroup);
			}
			groups.get(value).add(index);
		}
		Integer nM = groupsOfSizeAtLeastTwo.size();
		
		// Exclude fringe cases
		if (nM == 0) {
			// There is no group that could be split
			return Double.NEGATIVE_INFINITY;
		}
		
		// Generate the SPLIT
		Integer newIndex = parametersInput.get().size();
		// TODO: Some cleanup so this won't lead to a memory explosion 

		Integer groupToBeSplit = groupsOfSizeAtLeastTwo.get(Randomizer.nextInt(nM));
		
		ArrayList<Integer> oldGroup = groups.get(groupToBeSplit);
		Integer firstInNewGroup = Randomizer.nextInt(oldGroup.size()-1)+1;
		
		// Moving an entry from one group to another means changing the corresponding value in groupings.
		// At least one value has to move, so remove it from the oldGroup to not hit it twice.
		groupingsInput.get(this).get(oldGroup.remove((int) firstInNewGroup)).setValue(newIndex);
		// If we do not convert to (int), it tries to remove an element of that value,
		// as per the alternative definition of remove(Object O)
		
		Integer newGroupSize = 1;
		Integer oldGroupSize = oldGroup.size();
		
		// Go through the old list from the end, and either move or keep entries.
		for (int index = oldGroup.size()-1; index > 0; --index) {
			if (Randomizer.nextBoolean()) {
				groupingsInput.get(this).get(oldGroup.remove((int) index)).setValue(newIndex);
				++newGroupSize;
				--oldGroupSize;
			}
		}
		// Change the parameter value those entries now refer to, to reflect the old value.
		// TODO: Follow the Pagel & Meade paper, doing one of:
		//  * Calling an appropriate Operator
		//  * Implementing the change here
		//  * Making sure that it works without that step
		parametersInput.get(this).add(new RealParameter());
		parametersInput.get(this).get(newIndex).setValue(
				parametersInput.get(this).get(groupToBeSplit).getValue());
				
		// If only a split can happen, it has probability 1.
		double logSplitProbability = 0;
		// If splitting and merging can both happen, the split probability is 1/2.
		if (groups.size() != 1) {
			logSplitProbability = Math.log(0.5);
		}
		
		// If, after this, only a merge can happen, that merge has probability 1.
		double logMergeProbability = 0;
		// If splitting and merging will both be options, the merge probability is 1/2.
		// This is not the case if we split the last group of size at least two into two
		// single groups of size one.
		if (groupsOfSizeAtLeastTwo.size() != 1 || newGroupSize != 1 || oldGroupSize != 1) {
			logMergeProbability = Math.log(0.5);
		}
		
		// The proposal ratio for for a split move is
		// [ P_m(M') 1/(k nCr 2) ]/[ P_s(M) 1/N(M) 1/(2^(n_i+n_j-1)-1) 1/(q (n_i+n_j)) ]
		
		return logMergeProbability - Binomial.logChoose(k, 2)
				- logSplitProbability + Math.log(groupsOfSizeAtLeastTwo.size())
				+ Math.log(Math.pow(2, newGroupSize+oldGroupSize-1)-1)
				+ Math.log(
						// TODO: Understand how the rate plays a role here
						parametersInput.get().get(groupToBeSplit).getValue() *
						(newGroupSize+oldGroupSize));
	}
}
