# Practical 1: Rigid-Body Simulation

##Handout date: 20/Feb/2018.

##Deadline: 13/Mar/2018 9:00.

The first practical is about simulating a the motion of rigid bodies, and their interactions through collision. The objectives of the practical are:

1. Implement rigid velocity, position, and orientation integration in a discrete-time iteration.
 <br />
2. Implement interpenetration resolution.
 <br />
2. Implement impulse-based collision resolution.
 <br />
3. Extend the framework with some chosen effects.  

This is the repository for the skeleton on which you will build your first exercise. Using CMake allows you to work and submit your code in all platforms. The entire environment is in C++, but most of the "nastier" coding parts have been simplified; for the most part, you only code the mathemtical-physical parts.

![alt text](snap.png)
on as substitute to **one** in the list above, but it needs approval on the Lecturer's behalf **beforehand**.


##Installation

The skeleton uses the following dependencies: [libigl](http://libigl.github.io/libigl/), and consequently [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), for the representation and viewing of geometry, and [libccd](https://github.com/danfis/libccd) for collision detection. libigl viewer is using [imGUI](https://github.com/ocornut/imgui) for the menu. Everything is bundled as either submodules, or just incorporated code within the environment, and you do not have to take care of any installation details. To get the library, use:

```bash
git clone --recursive https://github.com/avaxman/INFOMGP-Practical1.git
```

to compile the environment, go into the `practical1` folder and enter in a terminal (macOS/Linux):

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

In windows, you need to use [cmake-gui](https://cmake.org/runningcmake/). Pressing twice ``configure`` and then ``generate`` will generate a Visual Studio solution in which you can work. The active soution should be ``practical1_bin``. *Note*: it only seems to work in 64-bit mode. 32-bit mode might give alignment errors.

##Using the dependencies

You do not need to utilize any dependency on your own, or install anything other than the above. For the most part, the dependencies are parts of code that are background, or collision detection code, which is not a direct part of the practical. The exception is [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for the representation and manipulation of vectors and matrices, but it is a quite a shallow learning curve. It is even possible to learn most necessary aspects from looking at the existing code. However, it is advised to go through the "getting started" section on the Eigen website (reading up to and including "Dense matrix and array manipulation" should be enough). Adding options in [imGUI](https://github.com/ocornut/imgui) is also required for some extensions; for that, their page gives many examples, and tutorial 106_ViewerMenu in [libigl](http://libigl.github.io/libigl/) also provides one specific o the GUI of this exercise.











