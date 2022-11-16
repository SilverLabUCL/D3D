package ucl.silver.d3d.core;

/**
 * <p>
 * Title: D3D</p>
 * <p>
 * Description: 3D Reaction-Diffusion Simulator</p>
 * <p>
 * Copyright: Copyright (c) 2022</p>
 * <p>
 * Company: The Silver Lab at University College London</p>
 *
 * @author Jason Rothman
 * @version 2.1
 */
public class D3Derror extends ParamVector {

    public String object = "";
    public String function = "";
    public String parameter = "";
    public String error = "";

    @Override
    public boolean canEdit(String name) {
        if (name.equalsIgnoreCase("object")) {
            return false;
        }
        if (name.equalsIgnoreCase("function")) {
            return false;
        }
        if (name.equalsIgnoreCase("parameter")) {
            return false;
        }
        if (name.equalsIgnoreCase("error")) {
            return false;
        }
        return true;
    }

    public D3Derror(String OBJECT, String FUNCTION, String PARAMETER, String ERROR) {
        super(null);
        object = OBJECT;
        function = FUNCTION;
        parameter = PARAMETER;
        error = ERROR;
        createVector(true);
    }

}
