package ucl.silver.d3d.gui;

import javax.swing.*;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeEvent;
import java.text.*;

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
public class UserInputField implements PropertyChangeListener {

    public String vName;
    public JLabel iLabel;
    public UserInputTextField iField;
    public NumberFormat iFormat;
    public UserInputAction action;

    public UserInputField(String n, String label, Double d, int col, UserInputAction ua) {

        action = ua;
        vName = n;
        iLabel = new JLabel(label);

        // create text field
        setUpFormats();
        iField = new UserInputTextField(iFormat);
        iField.setValue(d);
        iField.setColumns(col);
        iField.addPropertyChangeListener("value", this);
        iField.selectAll();

        // tell accessibility tools about label/textfield pairs
        iLabel.setLabelFor(iField);

    }

    public Number getUserValue() {
        return (Number) iField.getValue();
    }

    private void setUpFormats() {
        iFormat = NumberFormat.getNumberInstance();
    }

    /** Called when a field's "value" property changes. */
    public void propertyChange(PropertyChangeEvent e) {
        Object source = e.getSource();
        Number n = (Number) iField.getValue();
        if (source != iField) {
            return;
        }
        action.actionField(e, vName, n);
    }
}

