package ucl.silver.d3d.gui;

import javax.swing.JFormattedTextField;
import java.awt.event.FocusEvent;
import java.text.NumberFormat;

/**
 * <p>
 * Title: D3D</p>
 * <p>
 * Description: 3D Reaction-Diffusion Simulator</p>
 * <p>
 * Copyright: Copyright (c) 2018</p>
 * <p>
 * Company: The Silver Lab at University College London</p>
 *
 * @author Jason Rothman
 * @version 1.0
 */
public class UserInputTextField extends JFormattedTextField {

  UserInputTextField(NumberFormat nf) {
    super(nf);
  }

  @Override
  protected void processFocusEvent(FocusEvent e)
  {
    super.processFocusEvent(e);
    if (e.getID() == FocusEvent.FOCUS_GAINED)
    {
      this.selectAll();
    }
    else
    {
      this.select(0, 0);
    }
  }
}
