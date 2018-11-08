package ucl.silver.d3d.gui;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;
import java.util.*;
import java.text.DecimalFormat;
import ucl.silver.d3d.core.*;

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
public class PanelParamsTable extends JTable {

    Vector<ParamObject> dataVector = new Vector<ParamObject>(2, 2);

    Color green = new Color(153, 153, 102);
    Color grey = new Color(230, 230, 230);
    Color red = new Color(204, 0, 0);

    int cwidth0 = 200;
    int cwidth1 = 200;
    int cwidth2 = 50;
    int cwidth3 = 50;
    int cwidth4 = 50;
    int cwidth5 = 200;

    DecimalFormat dFormat = new DecimalFormat("0.0####");
    DecimalFormat dFormatSmall = new DecimalFormat("0.000E0");

    ParamVector paramVectorSelect = null;

    public PanelParamsTable() {

        super();

        setModel(new SimpleTableModel());

        TableColumn col_0 = getColumnModel().getColumn(0);
        col_0.setPreferredWidth(cwidth0);

        TableColumn col_1 = getColumnModel().getColumn(1);
        col_1.setPreferredWidth(cwidth1);

        TableColumn col_2 = getColumnModel().getColumn(2);
        col_2.setPreferredWidth(cwidth2);

        TableColumn col_3 = getColumnModel().getColumn(3);
        col_3.setPreferredWidth(cwidth3);

        TableColumn col_4 = getColumnModel().getColumn(4);
        col_4.setPreferredWidth(cwidth4);

        TableColumn col_5 = getColumnModel().getColumn(5);
        col_5.setPreferredWidth(cwidth5);

        //setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        //setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);

        //setIntercellSpacing(new Dimension(10, 1));

        //setMinimumSize(new Dimension(width(), 2000));
        //setPreferredSize(new Dimension(width(), 2000));

    }

    public int width(){
        return cwidth0 + cwidth1 + cwidth2 + cwidth3 + cwidth4 + cwidth5;
    }

    @Override
    public Component prepareRenderer(TableCellRenderer renderer, int Index_row, int Index_col) {

        Component comp = super.prepareRenderer(renderer, Index_row, Index_col);

        String parameter = getModel().getValueAt(Index_row, 0).toString();
        String valueStr = getModel().getValueAt(Index_row, 1).toString();
        String units = getModel().getValueAt(Index_row, 2).toString();
        String type = getModel().getValueAt(Index_row, 3).toString();

        Color color;

        ParamObject data;

        boolean canEdit = true;

        if (dataVector == null) {
            return comp;
        }

        data = (ParamObject) (dataVector.elementAt(Index_row));

        if (data.paramVector != null) {
            canEdit = data.paramVector.canEdit(data.getName());
        }

        if (data.ERROR != null) {
            comp.setBackground(red);
            comp.setForeground(Color.WHITE);
        } else if (type.matches("Class")) {
            comp.setBackground(green);
            comp.setForeground(Color.WHITE);
        } else if (units.matches("RGB")) {
            color = ColorD3D.string2color(valueStr);
            comp.setBackground(color);
            color = ColorD3D.getColorForeground(valueStr);
            comp.setForeground(color);
        } else if (type.equalsIgnoreCase("")) {
            comp.setBackground(Color.white);
            comp.setForeground(Color.BLACK);
        } else if (!canEdit) {
            comp.setBackground(grey);
            comp.setForeground(Color.BLACK);
        } else {
            comp.setBackground(Color.white);
            comp.setForeground(Color.BLACK);
        }

        return comp;

    }

    // function to create parameter vector from param table
    public void setFields() {

        ParamObject[] v = null;

        if (paramVectorSelect != null) {
            paramVectorSelect.updateVectors();
            v = paramVectorSelect.getVector();
        }

        dataVector.removeAllElements();

        if (v != null) {

            for (ParamObject o : v) {

                if (o.getType().equalsIgnoreCase("EndClass")) {
                    continue;
                }

                dataVector.addElement(o);

            }

        }

        revalidate();
        repaint();

    }

    // inner class
    public class SimpleTableModel extends AbstractTableModel {

        public String[] colNames = {"Parameter", "Value", "Units", "Type", "Error", "Definition"}; // column name
        public Class[] colClass = {String.class, String.class, String.class, String.class, String.class, String.class}; // column type

        @Override
        public int getColumnCount() {
            return colNames.length;
        }

        @Override
        public int getRowCount() {
            return dataVector.size();
        }

        /**
         * setValueAt
         * This function updates the data in the TableModel
         * depending upon the change in the JTable
         */
        @Override
        public void setValueAt(Object value, int row, int col) {

            String newColorStr, fileName;
            boolean success = false;

            if (dataVector == null) {
                return;
            }

            ParamObject data = (ParamObject) (dataVector.elementAt(row));

            if ((col == 1) && (data != null) && (paramVectorSelect != null)) {

                if (data.getUnits().equalsIgnoreCase("RGB")) {

                    newColorStr = ColorD3D.promptColorOrScale(data.getValue());

                    value = newColorStr;

                    if (newColorStr != null) {

                        success = paramVectorSelect.set(data, newColorStr);

                        if (!success) {
                            if (data.paramVector != null) {
                                success = data.paramVector.set(data, newColorStr);
                            }
                        }

                    }
                    
                } else if (data.getUnits().equalsIgnoreCase("DIR")) {

                    fileName = Master.getFileName(value.toString());
                    
                    success = paramVectorSelect.set(data, fileName);

                    if (!success) {
                        if (data.paramVector != null) {
                            success = data.paramVector.set(data, fileName);
                        }
                    }

                } else {

                    success = paramVectorSelect.set(data, value.toString());

                    if (!success) {
                        if (data.paramVector != null) {
                            success = data.paramVector.set(data, value.toString());
                        }
                    }

                }

                if (success) {
                    Master.log(paramVectorSelect.getClass().getSimpleName() + "." + data.getName() + " = " + data.getValue());
                }

            }

            revalidate();
            repaint();

        }

        @Override
        public String getColumnName(int col) {
            return colNames[col];
        }

        @Override
        public Class getColumnClass(int col) {
            return colClass[col];
        }

        /**
         * getValueAt
         * This function updates the JTable depending upon the
         * data in the TableModel
         */
        @Override
        public Object getValueAt(int row, int col) {

            ParamObject data = (ParamObject) (dataVector.elementAt(row));

            String type, value;

            double d;

            switch (col) {

                case 0:
                    return data.getName();

                case 1:

                    type = data.getType();
                    value = data.getValue();

                    if (type.equalsIgnoreCase("double")) {

                        d = Double.parseDouble(value);

                        if (d == 0) {
                            return value;
                        } else if (Math.abs(d) < 0.0001) {
                            return dFormatSmall.format(d);
                        } else {
                            return dFormat.format(d);
                        }

                    }

                    return value;

                case 2:
                    return data.getUnits();

                case 3:
                    return data.getType();

                case 4:
                    return data.ERROR;

                case 5:
                    return data.getHelp();

            }

            return new String();

        }

        @Override
        public boolean isCellEditable(int row, int col) {

            if (dataVector == null) {
                return false;
            }

            ParamObject data = (ParamObject) (dataVector.elementAt(row));

            if (data.paramVector == null) {
                return false;
            }

            boolean canEdit = data.paramVector.canEdit(data.getName());

            if (col == 1) {
                return canEdit;
            }

            return false;

        }
    }
}
