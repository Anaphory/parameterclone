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

package parameterclone.selector;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

@Description("A calculation node that propagates parameters from a vector of parameters")
// Needed for the reversible-jump Markov chain described in
@Citation("Pagel, M., Meade, A., 2006. Bayesian Analysis of Correlated Evolution of Discrete Characters by Reversible-Jump Markov Chain Monte Carlo. The American Naturalist 167, 808--825. doi:10.1086/503444")
public class Selector extends CalculationNode implements Loggable, Function {
	// Input objects
	final public Input<IntegerParameter> entryInput = new Input<IntegerParameter>(
			"entry",
			"The index of the parameter vector that this object propagates",
			// TODO: This should be made optional, defaulting to
			// 0..parameters.getDimension()
			Validate.REQUIRED);
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			new RealParameter(), Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices", new IntegerParameter(),
			Validate.REQUIRED);
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each parameter",
			Validate.REQUIRED);

	// Member objects
	protected IntegerParameter entries;
	protected Integer maxIndex;

	@Override
	public void initAndValidate() throws Exception {
		maxIndex = parametersInput.get().getDimension();
		if (entryInput == null) {
			// fill it with 0..maxIndex-1
		} else {
			entries = entryInput.get();
			for (Integer entry : entries.getValues()) {
				if (entry > groupingsInput.get().getDimension()) {
					throw new Exception(
							"entries must be valid index of groupings");
				}
			}
		}
		// Array-like RealParameters do not implement java.lang.iterable, so we
		// must do the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new Exception(
						"All entries in groupings must be valid indices of parameters");
			}
		}
		// value = parametersInput[groupingsInput[entry]]
	}

	/**
	 * Function interface implementation follows *
	 */

	@Override
	public int getDimension() {
		return entries.getDimension();
	}

	@Override
	public double getArrayValue() {
		return getArrayValue(0);
	}

	@Override
	public double getArrayValue(int iDim) {
		int index = groupingsInput.get().getNativeValue(
				entries.getNativeValue(iDim));
		return parametersInput.get().getValue(index);
	}

	/**
	 * Loggable interface implementation follows *
	 */

	@Override
	public void init(final PrintStream out) throws Exception {
		for (int i = 0; i < getDimension(); ++i) {
			out.print(getID() + "" + i + "\t");
		}
	}

	@Override
	public void log(final int nSample, final PrintStream out) {
		for (int i = 0; i < getDimension(); ++i) {
			out.print(getArrayValue(i) + "\t");
		}
	}

	@Override
	public void close(final PrintStream out) {
		// nothing to do
	}

}
