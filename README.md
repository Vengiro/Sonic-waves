# MATLAB DSP Software for Wireless Image Transmission

## Overview
This project implements a complete Digital Signal Processing (DSP) system in MATLAB for wireless image transmission and reception. It simulates the processes involved in encoding an image into a signal, transmitting it via a 200 MHz carrier, and reconstructing the original image from the received signal. The system adheres to time and power constraints while employing advanced modulation and signal processing techniques to ensure reliable communication.

---

## Features
1. **Image-to-Signal Transformation**:
   - Converts an image into binary data.
   - Encodes binary data into signals using various modulation schemes.
   - Ensures compliance with specified time and power constraints.

2. **Signal Transmission**:
   - Transforms the processed signal into a complex baseband signal.
   - Uploads the signal to a radio transmitter operating at 200 MHz.
   - Sends the signal to a paired radio receiver.

3. **Signal Reception and Image Reconstruction**:
   - Processes the raw received signal to reconstruct the original image.
   - Implements methods for:
     - Frame recovery
     - Timing recovery
     - Equalization

4. **Modulation Techniques**:
   - Basic Modulations: BPSK, 16-QAM.
   - Advanced Modulation: Trellis Coded Modulation (TCM) with a 16-PSK constellation.

5. **Performance Enhancements**:
   - Robust against noise and channel impairments.
   - Includes advanced recovery techniques for precise signal decoding.

---

## System Components

### Transmitter
1. **Image Processing**:
   - Reads the input image and converts it into binary data.
   - Optionally compresses or resizes the image for optimized transmission.

2. **Signal Encoding**:
   - Maps binary data to modulation symbols using selected modulation schemes (BPSK, 16-QAM, or 16-PSK).
   - Applies pulse shaping to meet time and power constraints.

3. **Radio Upload**:
   - Converts the baseband signal to an appropriate format.
   - Transmits the signal at a 200 MHz carrier frequency.

### Receiver
1. **Signal Acquisition**:
   - Receives the transmitted signal via a paired radio.
   - Converts the RF signal to a complex baseband format for processing.

2. **Signal Decoding**:
   - Performs frame and timing recovery to extract data frames.
   - Equalizes the channel to mitigate distortions.
   - Decodes modulation symbols back to binary data.

3. **Image Reconstruction**:
   - Reconstructs the original image from binary data.
   - Displays the final image for verification.

---

## Advanced Techniques
1. **Trellis Coded Modulation (TCM)**:
   - Enhances error performance using convolutional encoding.
   - Utilizes a 16-PSK constellation for higher spectral efficiency.

2. **Equalization**:
   - Employs adaptive algorithms to compensate for channel impairments.

3. **Timing Recovery**:
   - Implements robust methods to align transmitted and received signals.

---

## Requirements
- **MATLAB**: R2021b or later.
- **Hardware**: A pair of radios capable of 200 MHz transmission and reception.
- **Image Input**: BMP images of any size.

---

## Usage
1. **Setup**:
   - Install MATLAB and necessary toolboxes.
   - Configure the radio hardware for 200 MHz operation.

2. **Running the Code**:
   - Modify the input parameters in the MATLAB script (e.g., image file, modulation scheme).
   - Execute the transmitter script to upload the signal to the radio.
   - Run the receiver script to process the received signal.

3. **Results**:
   - View the reconstructed image in MATLAB.
   - Analyze performance metrics (e.g., Bit Error Rate, Signal-to-Noise Ratio).

---

## Customization
- **Modulation Schemes**: Switch between BPSK, 16-QAM, and 16-PSK in the script.
- **Advanced Techniques**: Enable or disable TCM and recovery algorithms based on system requirements.
- **Performance Tuning**: Adjust parameters such as symbol rate, power levels, and filter designs for optimal results.

---

## Troubleshooting
- Ensure proper synchronization between the transmitter and receiver radios.
- Verify that the radio hardware meets the required specifications.
- Debug frame and timing recovery stages if the reconstructed image is distorted.

---

## Future Improvements
- Extend to higher-order modulation schemes (e.g., 64-QAM).
- Implement additional error correction methods.
- Enhance system robustness under varying channel conditions.

---

## Acknowledgments
This project combines theoretical DSP knowledge and practical wireless communication techniques to achieve efficient image transmission. It serves as a foundation for exploring advanced wireless system designs.

