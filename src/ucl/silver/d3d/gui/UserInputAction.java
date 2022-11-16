package ucl.silver.d3d.gui;

import java.beans.PropertyChangeEvent;
import java.awt.event.ActionEvent;

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
public interface UserInputAction {
  public void actionField(PropertyChangeEvent e, String name, Number n);
  public void actionButton(ActionEvent e);
}
