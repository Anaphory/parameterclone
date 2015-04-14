/**
 * 
 */
package parameterclone.selectedprior;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.Prior;

/**
 * @author gereon
 *
 */
public class SelectedPrior extends Prior {
	public Input<BooleanParameter> active = new Input<BooleanParameter>("active", "array denoting whether dimensions of x are active", Validate.REQUIRED);

	ParametricDistribution dist;

	@Override
    public void initAndValidate() throws Exception {
		if (!(m_x.get() instanceof RealParameter)) {
			throw new Exception(
					"SelectedPrior only takes RealParameter arguments, because it is intended to work with Selectors");
		}
    	if (active.get().getDimension() != m_x.get().getDimension()) {
    		throw new Exception(
    				"Dimension of active markers does not correspond to dimension of x");
    	}
        dist = distInput.get();
        calculateLogP();
    }

    @Override
    public double calculateLogP() throws Exception {
        Function x = m_x.get();
        
        // cut x down to include only the active entries of m_x
        Integer sum = 0;
        for (Boolean i : active.get().getValues()) {
            sum += i ? 1 : 0;
        }
        
        // essentially, cut_x = m_x[active]
        RealParameter cut_x = new RealParameter();
        cut_x.setDimension(sum);
        int j = sum - 1;
        for (int i = x.getDimension() - 1; i>=0; --i) {
        	if (active.get().getValue(i)) {
            	cut_x.setValue(j, x.getArrayValue(i));
            	--j;
			}
		}
        
		// test that parameter is inside its bounds
		double l = 0.0;
		double h = 0.0;
		l = ((RealParameter) x).getLower();
		h = ((RealParameter) x).getUpper();

		for (int i = 0; i < cut_x.getDimension(); i++) {
			double value = cut_x.getArrayValue(i);
			if (value < l || value > h) {
				logP = Double.NEGATIVE_INFINITY;
				return Double.NEGATIVE_INFINITY;
			}
		}

		logP = dist.calcLogP(cut_x);
        return logP;
    }

}
