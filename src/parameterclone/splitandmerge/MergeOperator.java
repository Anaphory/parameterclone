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
import java.util.HashSet;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.Binomial;
import beast.util.Randomizer;

@Description("Randomly merge two groups of parameters")
@Citation("Pagel, M., Meade, A., 2006."
		+ " Bayesian Analysis of Correlated Evolution of Discrete Characters"
		+ " by Reversible-Jump Markov Chain Monte Carlo."
		+ " The American Naturalist 167, 808--825. doi:10.1086/503444")
public class MergeOperator extends Operator {
	// Inputs that are changed by the operator
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices", Validate.REQUIRED);
	public Input<BooleanParameter> parameterIsUsed = new Input<BooleanParameter>(
			"parameterIsUsed",
			"stores whether the corresponding element of parameters is active",
			Validate.REQUIRED);

	Integer k;
	Integer maxIndex;

	@Override
	public void initAndValidate() throws Exception {
		k = groupingsInput.get().getDimension();
		maxIndex = parametersInput.get().getDimension();
		// Array-like RealParameters do not implement java.lang.iterable, so we
		// must do the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new Exception(
						"All entries in groupings must be valid indices of parameters");
			}
		}
		if (parameterIsUsed.get().getDimension() != parametersInput.get()
				.getDimension()) {
			throw new Exception(
					"parameterIsUsed must correspond to parameters in dimension");
		}
	}

	/**
	 * Change the parameter and return the log of the Hastings ratio: Merge two
	 * groups of joined parameters, averaging the parameter.
	 */
	@Override
	public double proposal() {

		try {

			// Find the composition of groups

			// For each grouping value, the corresponding indices of groupings
			HashMap<Integer, HashSet<Integer>> groups = new HashMap<>();
			ArrayList<Integer> groupIndices = new ArrayList<>(); // Keys of
																	// groups
			HashSet<Integer> indicesOccuringAtLeastTwice = new HashSet<>(); // large
																			// groups

			for (Integer index = groupingsInput.get().getDimension() - 1; index >= 0; --index) {
				Integer value = groupingsInput.get().getNativeValue(index);
				if (groups.get(value) != null) {
					if (groups.containsKey(value)) {
						indicesOccuringAtLeastTwice.add(value);
					}
				} else {
					groupIndices.add(value);
					HashSet<Integer> newGroup = new HashSet<>();
					groups.put(value, newGroup);
				}
				groups.get(value).add(index);
			}

			Integer nGroups = groupIndices.size();

			Integer rawMergeIndex = Randomizer.nextInt(nGroups);
			Integer rawRemoveIndex = Randomizer.nextInt(nGroups - 1);
			if (rawRemoveIndex >= rawMergeIndex) {
				++rawRemoveIndex;
			}
			Integer mergeIndex = groupIndices.get(Randomizer.nextInt(nGroups));
			Integer removeIndex = groupIndices.get(Randomizer.nextInt(nGroups));

			HashSet<Integer> mergeGroup = groups.get(mergeIndex);
			HashSet<Integer> removeGroup = groups.get(removeIndex);

			Integer mergeGroupSize = mergeGroup.size();
			Integer removeGroupSize = removeGroup.size();

			double newValue = (parametersInput.get(this).getValue(mergeIndex)
					* mergeGroupSize + parametersInput.get(this).getValue(
					removeIndex)
					* removeGroupSize)
					/ (mergeGroupSize + removeGroupSize);

			// Generate the MERGE

			parametersInput.get(this).setValue(mergeIndex, newValue);
			for (Integer toBeMerged : removeGroup) {
				// groupings[toBeMerged] = mergeIndex
				groupingsInput.get(this).setValue(toBeMerged, mergeIndex);
			}
			parameterIsUsed.get(this).setValue(removeIndex, false);

			// Now we calculate the Hastings ratio.

			Integer groupsOfSizeAtLeastTwo = indicesOccuringAtLeastTwice.size();

			// If only a merge can happen, it has probability 1.
			double logMergeProbability = 0;
			// If splitting and merging can both happen, the merge probability
			// is 1/2.
			// Splitting can happen if
			if (groupsOfSizeAtLeastTwo != 0) {
				logMergeProbability = Math.log(0.5);
			}

			// If we merged two groups of size one, we gain a group of size at
			// least two.
			if (mergeGroupSize == 1 && removeGroupSize == 1) {
				++groupsOfSizeAtLeastTwo;
			}
			// If we merged two groups of size at least two, we lose one in
			// number.
			if (mergeGroupSize >= 2 && removeGroupSize >= 2) {
				--groupsOfSizeAtLeastTwo;
			}

			// If, after this, only a split can happen, that split has
			// probability 1.
			double logSplitProbability = 0;
			// If splitting and merging will both be options, the split
			// probability is 1/2.
			// This is not the case if we merged the last two groups.
			if (nGroups != 2) {
				logSplitProbability = Math.log(0.5);
			}

			// The proposal ratio for for a merge move is
			// [ P_s(M') 1/N(M') 1/(2^(n'_i+n'_j-1)-1) 1/(q' (n'_i+n'_j)) ]/[
			// P_m(M) 1/(k nCr 2) ]

			return logSplitProbability
					- Math.log(groupsOfSizeAtLeastTwo)
					- Math.log(Math
							.pow(2, mergeGroupSize + removeGroupSize - 1) - 1)
					- Math.log(
					// TODO: Understand how the rate plays a role here
					parametersInput.get().getValue(mergeIndex)
							* (mergeGroupSize + removeGroupSize))
					- logMergeProbability + Binomial.logChoose(k, 2);

		} catch (Exception e) {
			// whatever went wrong, we want to abort this operation...
			return Double.NEGATIVE_INFINITY;
		}

	}
}
