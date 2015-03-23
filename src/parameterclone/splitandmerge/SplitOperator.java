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
import java.util.Vector;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.math.Binomial;
import beast.util.Randomizer;

@Description("Randomly split a group of parameters in two")
@Citation("Pagel, M., Meade, A., 2006."
		+" Bayesian Analysis of Correlated Evolution of Discrete Characters"
		+" by Reversible-Jump Markov Chain Monte Carlo."
		+" The American Naturalist 167, 808--825. doi:10.1086/503444")
public class SplitOperator<T> extends Operator {
	// Inputs that are changed by the operator
	public Vector<Input<Integer>> groupings =
			new Vector<Input<Integer>>();//groupings -- A vector assigning parameter array indices to entries
	public Vector<Input<T>> parameters =
			new Vector<Input<T>>();//parameters -- A vector of connected input parameters

	Integer k;

	@Override
	public void initAndValidate() throws Exception {
		k = parameters.size();
		for (Input<Integer> value : groupings) {
			if (value.get() >= k) {
				throw new Exception("All entries in groupings must be valid indices of parameters.");
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
		for (Integer index = 0; index <= groupings.size(); ++index) {
			Integer value = groupings.get(index).get(this);
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
		if (nM >= k) {
			// There is no parameter value left to split to
			// TODO: Alternatively, we can invent a new parameter value. Would that make more sense?
			return Double.NEGATIVE_INFINITY;
		}
		
		// Generate the SPLIT
		Integer newIndex = null;
		for (Integer i = 0; i < groupings.size(); ++i) {
			if (!groupings.contains(i)) {
				newIndex = i;
				break;
			}
		}
		if (newIndex == null) {
			// This should not have happened, we should be able to name another group.
			return Double.NEGATIVE_INFINITY;
		}

		Integer groupToBeSplit = groupsOfSizeAtLeastTwo.get(Randomizer.nextInt(nM));
		
		ArrayList<Integer> oldGroup = groups.get(groupToBeSplit);
		Integer firstInNewGroup = Randomizer.nextInt(oldGroup.size()-1)+1;
		
		// Moving an entry from one group to another means changing the corresponding value in groupings.
		// At least one value has to move, so remove it from the oldGroup to not hit it twice.
		groupings.get(oldGroup.remove((int) firstInNewGroup)).setValue(newIndex, this);
		// If we do not convert to (int), it tries to remove an element of that value,
		// as per the alternative definition of remove(Object O)
		
		Integer newGroupSize = 1;
		Integer oldGroupSize = oldGroup.size();
		
		// Go through the old list from the end, and either move or keep entries.
		for (int index = oldGroup.size()-1; index > 0; --index) {
			if (Randomizer.nextBoolean()) {
				groupings.get(index).setValue(newIndex, this);
				++newGroupSize;
				--oldGroupSize;
			}
		}
		// Change the parameter value those entries now refer to, to reflect the old value.
		// TODO: Follow the Pagel&Meade paper, doing one of:
		//  * Calling an UpDownOperator
		//  * Implementing the change here
		//  * Making sure that it works without that step
		parameters.get(newIndex).setValue(parameters.get(groupToBeSplit).get(this), this);
		
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
						// parameters.get(groupToBeSplit).get(this) *
						(newGroupSize+oldGroupSize));
	}
}
