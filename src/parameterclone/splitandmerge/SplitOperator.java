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

@Description("Randomly split a group of parameters in two")
@Citation("Pagel, M., Meade, A., 2006."
		+ " Bayesian Analysis of Correlated Evolution of Discrete Characters"
		+ " by Reversible-Jump Markov Chain Monte Carlo."
		+ " The American Naturalist 167, 808--825. doi:10.1086/503444")
public class SplitOperator extends Operator {
	// Inputs that are changed by the operator
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			new RealParameter(), Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices", new IntegerParameter(),
			Validate.REQUIRED);
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each parameter",
			(IntegerParameter) null);

	Integer maxIndex;

	@Override
	public void initAndValidate() throws Exception {
		maxIndex = parametersInput.get().getDimension();
		// RealParameter does not implement java.lang.iterable, so we must do
		// the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new Exception(
						"All entries in groupings must be valid indices of parameters");
			}
		}
	}

	/**
	 * Change the parameter and return the log of the Hastings ratio. Split a
	 * class of joined parameters in two.
	 */
	@Override
	public double proposal() {
		// Find the composition of groups, in particular which ones can be
		// split.

		// If only a split can happen, it has probability 1.
		// If splitting and merging can both happen, the split probability
		// is 1/2.
		Integer groups = 0;
		Integer nM = 0;
		for (int i = parametersInput.get(this).getDimension() - 1; i >= 0; --i) {
			Integer size = sizesInput.get(this).getNativeValue(i);
			if (size > 0) {
				++groups;
				if (size >= 2) {
					nM += 1;
				}
			}
		}
		double logSplitProbability = groups > 1 ? Math.log(0.5) : 0.;

		// Exclude fringe cases
		if (nM == 0) {
			// There is no group that could be split
			// System.out.printf("Split -- No group to be split\n");
			return Double.NEGATIVE_INFINITY;
		}

		// Generate the parameter index for the new group and find the index of
		// the group to be split
		Integer rawGroupToBeSplit = Randomizer.nextInt(nM);
		Integer groupToBeSplit = 0;
		Integer newIndex = null;
		Integer steps = sizesInput.get(this).getDimension();
		Integer k = 0;
		for (int index = 0; index < steps; ++index) {
			// Check the available indices from low to high:
			if (sizesInput.get(this).getValue(index) == 0) {
				// use the first available index encountered
				if (newIndex == null) {
					newIndex = index;
				}
			} else {
				++k;
			}
			if (sizesInput.get(this).getValue(index) < 2
					& rawGroupToBeSplit >= 0) {
				++groupToBeSplit;
			} else {
				++groupToBeSplit;
				--rawGroupToBeSplit;
			}
		}
		groupToBeSplit = groupToBeSplit + rawGroupToBeSplit;
		if (newIndex == null) {
			// System.out.printf("Split -- No newIndex\n");
			return Double.NEGATIVE_INFINITY;
		}

		// Generate the SPLIT

		ArrayList<Integer> oldGroup = new ArrayList<Integer>();
		for (int index = groupingsInput.get(this).getDimension() - 1; index >= 0; --index) {
			if (groupingsInput.get(this).getNativeValue(index) == groupToBeSplit) {
				oldGroup.add(index);
			}
		}
		if (oldGroup.size() <= 1) {
			// System.out.printf("Split -- Splitting a non-group\n");
			return Double.NEGATIVE_INFINITY;
		}

		// System.out.printf("Split %d into %d\n", groupToBeSplit, newIndex);

		// Moving an entry from one group to another means changing the
		// corresponding value in groupings.
		// I theory, there should be a way without rejections.
		// But for correct-before-efficient reasons, just generate any
		// partition, and reject the operator when the partition is trivial.

		Integer newGroupSize = 0;
		Integer oldGroupSize = oldGroup.size();

		// Go through the old list from the end, and either move or keep
		// entries. Note that index 0 is definitely staying in the old
		// group, so we only iterate while index > 0.
		for (int index = oldGroup.size() - 1; index > 0; --index) {
			if (Randomizer.nextBoolean()) {
				groupingsInput.get(this).setValue(oldGroup.remove((int) index),
						newIndex);
				++newGroupSize;
				--oldGroupSize;
			}
		}
		
		if (newGroupSize == 0 || oldGroupSize == 0) {
			// System.out.printf("Split -- Non-Split generated\n");
			return Double.NEGATIVE_INFINITY;
		}

		sizesInput.get(this).setValue(newIndex, newGroupSize);
		sizesInput.get(this).setValue(groupToBeSplit, oldGroupSize);

		Double logJacobian = -Math.log(newGroupSize + oldGroupSize)
				+ Math.log(newGroupSize) + Math.log(oldGroupSize);

		// Change the parameter value those entries now refer to, to reflect
		// the old value.
		// TODO: Follow the Pagel & Meade paper, doing one of:
		// * Calling an appropriate Operator
		// * Implementing the change here
		// * Making sure that it works without that step
		parametersInput.get(this).setValue(newIndex,
				parametersInput.get(this).getValue(groupToBeSplit));

		// In order to keep dimensions matched (cf. Green 1995, p. 716), there
		// needs to be a bijection between the pre-image and the image of this
		// operator and its inverse. This is mitigated by a random distortion of
		// the rates, keeping the sum of rates constant.
		// The proposal ration needs to take that into account.
		double rate = parametersInput.get(this).getValue(groupToBeSplit);
		double mu = Randomizer.uniform(-oldGroupSize * rate, newGroupSize
				* rate);
		parametersInput.get(this).setValue(groupToBeSplit,
				rate + mu / oldGroupSize);
		parametersInput.get(this).setValue(newIndex, rate - mu / newGroupSize);
		double bijectionDensity = Math
				.log(rate * (oldGroupSize + newGroupSize));

		// If, after this, only a merge can happen, that merge has
		// probability 1.
		// This is the case if we split the last group of size at least
		// two into two single groups of size one.
		// If splitting and merging will both be options, the merge
		// probability is 1/2.
		double logMergeProbability = (nM == 1 && newGroupSize == 1 && oldGroupSize == 1) ? 0
				: Math.log(0.5);

		// The proposal ratio for for a split move is
		// [ P_m(M') 1/(k' nCr 2) ]/[ P_s(M) 1/N(M) 1/(2^(n_i+n_j-1)-1) 1/(q
		// (n_i+n_j)) ]
		// NOTE: The reference states (k nCr 2), but that seems to be a typo. We
		// use k' = k+1 after a split.
		Double p = logMergeProbability - Binomial.logChoose(k + 1, 2)
				- logSplitProbability + Math.log(nM)
				+ Math.log(Math.pow(2, newGroupSize + oldGroupSize - 1) - 1)
				+ bijectionDensity + logJacobian;
		// System.out.printf("Split: %f\n", p);
		return p;
	}
}
