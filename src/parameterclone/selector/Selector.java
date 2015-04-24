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

@Description("A calculation node that propagates one parameter from a vector of parameters")
// Needed for the reversible-jump Markov chain described in
@Citation("Pagel, M., Meade, A., 2006. Bayesian Analysis of Correlated Evolution of Discrete Characters by Reversible-Jump Markov Chain Monte Carlo. The American Naturalist 167, 808--825. doi:10.1086/503444")
public class Selector extends CalculationNode implements Loggable, Function {
	// Input objects
	final public Input<Integer> entryInput = new Input<Integer>("entry",
			"The index of the parameter vector that this object propagates",
			Validate.REQUIRED);
	public Input<RealParameter> parametersInput = new Input<RealParameter>(
			"parameters",
			"individual parameters that the actual value is chosen from",
			new RealParameter(), Validate.REQUIRED);
	public Input<IntegerParameter> groupingsInput = new Input<IntegerParameter>(
			"groupings", "parameter selection indices",
			new IntegerParameter(), Validate.REQUIRED);
	public Input<IntegerParameter> sizesInput = new Input<IntegerParameter>(
			"sizes", "stores how many indices are pointing to each parameter",
			Validate.REQUIRED);

	// Member objects
	Integer entry;
	Integer maxIndex;
	Double value;
	Double storedValue;

	@Override
	public void initAndValidate() throws Exception {
		maxIndex = parametersInput.get().getDimension();
		entry = entryInput.get();
		if (entry > groupingsInput.get().getDimension()) {
			throw new Exception("entry must be valid index of groupings");
		}
		// Array-like RealParameters do not implement java.lang.iterable, so we must do the iteration by hand.
		for (int groupIndex = groupingsInput.get().getDimension() - 1; groupIndex >= 0; --groupIndex) {
			if (groupingsInput.get().getNativeValue(groupIndex) >= maxIndex) {
				throw new Exception(
						"All entries in groupings must be valid indices of parameters");
			}
		}
		// value = parametersInput[groupingsInput[entry]]
		Integer index = groupingsInput.get().getNativeValue(entry);
		value = parametersInput.get().getValue(index);
	}

	@Override
	public boolean requiresRecalculation() {
		// Essentially parameters[groupings[entry]], but wrapped as Inputs, Parameters and Lists
		value = parametersInput.get().getValue(groupingsInput.get().getNativeValue(entry));
		// Yes, this means we do the 'calculation' here in the check, but
		// it probably means avoiding other recalculations down the line.
		return storedValue != value;
	}

	@Override
	protected void store() {
		storedValue = value;
		super.store();
	}

	@Override
	protected void restore() {
		value = storedValue;
		storedValue = null;
		super.restore();
	}

	@Override
	protected void accept() {
		storedValue = null;
		super.accept();
	}

	/**
	 * Function interface implementation follows *
	 */

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue() {
		return value;
	}

	@Override
	public double getArrayValue(int iDim) {
		if (iDim == 0)
			return value;
		return 0;
	}

	/**
	 * Loggable interface implementation follows *
	 */

	@Override
	public void init(final PrintStream out) throws Exception {
		out.print(getID() + "\t");
	}

	@Override
	public void log(final int nSample, final PrintStream out) {
		out.print(value + "\t");
	}

	@Override
	public void close(final PrintStream out) {
		// nothing to do
	}

}
