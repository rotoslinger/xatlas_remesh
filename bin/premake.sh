#!/bin/bash

# Generate build files for macOS and Linux
echo "Generating build files..."

# Linux builds
premake5 gmake
premake5 --cc=clang gmake

# macOS builds
premake5 xcode4

echo "Build files generated successfully!"
echo ""
echo "To build on Linux:"
echo "  cd build/gmake"
echo "  make"
echo ""
echo "To build on macOS:"
echo "  open build/xcode4/xatlas.xcworkspace"
echo "  or use: xcodebuild -workspace build/xcode4/xatlas.xcworkspace -scheme viewer"
