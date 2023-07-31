# 5G Localization

This is a demonstration of 5G mmWave mMIMO localization.

This repository contains the functions and scripts necessary to test a geometry-based Mobile Station (MS) localization method in 5G networks.
 The method consists of three main stages:
1. **Channel modeling using**:
   - _mmWave_;
   - _mMIMO._
3. **Localization parameter estimation:**
   - _Time of Arrival (Toa);_
   - _2D Angle of Departure (AoD);_
   - _2D Angle of Arrival (AoA)._
4. **Localization:**
   - _Geometry-based algorithm for Line-of-Sight (LoS) condition;_
   - _Geometry-based algorithm for Non-Line-of-Sight (LoS) condition._

To perform the test, use the Matlab software and run the script **Main.m**. The default simulation is parameterized for LoS condition. However, it is possible to modify the parameters in the INPUT section of the **Main.m** script and conduct customized tests.
