package ucl.silver.d3d.core;

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
public class NativeC {

  // Loads the file NativeC.DLL at run-time
  static {
    try {
      System.loadLibrary("NativeC");
      System.out.println("Loaded library NativeC.");
    }
    catch (UnsatisfiedLinkError ule) {
      System.out.println("Cannot locate NativeC library.");
      //ule.printStackTrace();
    }
  }

  native public void displayHelloWorld();

  native public void integrateTest();

  //native public void integrateBasic(Project p, Compact c, boolean toFile);

  //native public void integrateUncage(Project p, Compact c, boolean toFile);

  // Initializes an array of elements 'count' times with value 'value'
  native public void initializeByteArray(byte[] array, int count, byte value);

  // Pass a 2-D array of bytes to be printed and modified
  native public void pass2DByteArray(byte[][] array);

  // Prints all members of this class by using callbacks to the method
  // 'toString' below
  native public String toStringWithPrint();

  // Prints a Java String in native code and returns a Java String
  native public String printLine(String text);

  // Prints w, x, y, and z using native code (no callback to this class)
  native public void printWXYZ();

  // Causes an exception, doesn't handle it in native code
  native public void causeException();

  // Causes an exception, handles it in native code
  native public void handleException();

  // Prints each element of the objArray using a callback to the object's
  // 'toString' method
  native public void printObjectArray(Object[] objArray, boolean Print);

  // Creates, allocates, and initializes an array or Rectangles in native
  // code and returns the array to Java
  native public java.awt.Rectangle[] returnRectArray(int size);

  // No arguments, no return value
  native public void VoidVoid();

  // For timing test of C++ calling back to Java
  native public void callbackVoid(int Count);

}
