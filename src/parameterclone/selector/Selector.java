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

import java.util.Vector;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;

@Description("A calculation node that propagates one parameter from a vector of parameters")
// Needed for the reversible-jump Markov chain described in
@Citation("Pagel, M., Meade, A., 2006. Bayesian Analysis of Correlated Evolution of Discrete Characters by Reversible-Jump Markov Chain Monte Carlo. The American Naturalist 167, 808--825. doi:10.1086/503444")
public class Selector<T> extends CalculationNode {
	// Input objects
	public Input<Integer> entryInput = new Input<Integer>(
			"entry","The index of the parameter vector that this object propagates",
			Validate.REQUIRED);
	public Vector<Input<T>> parameters =
			new Vector<Input<T>>();//parameters -- A vector of connected input parameters
	public Vector<Input<Integer>> groupings =
			new Vector<Input<Integer>>();//groupings -- A vector assigning each possible value of entry an index into parameters

	// Member objects
	Integer entry;
	
	@Override
	public void initAndValidate() throws Exception {
		entry = entryInput.get();
		
		if (entry >= groupings.size()) {
			throw new Exception("entry must be a valid index into groupings");
		}
		if (groupings.get(entry).get() >= parameters.size()) {
			throw new Exception("groupings values must be a valid index into parameters");
		}
	}
	
	public Input<T> get() {
		return parameters.get(groupings.get(entry).get());
	}
	
	@Override
	public boolean requiresRecalculation() {
		if (groupings.get(entry).isDirty()) {
			return true;
		}
		if (parameters.get(groupings.get(entry).get()).isDirty()) {
			return true;
		}
		return false;
	}

}
