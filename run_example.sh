#!/bin/bash

# xatlas Example Runner Script
# Usage: ./run_example.sh [model_name] [options]
# Examples:
#   ./run_example.sh cube
#   ./run_example.sh gazebo -verbose
#   ./run_example.sh Gargoyle_quadriflow

# Default values
EXAMPLE_BIN="build/gmake_clang/bin/x86_64/Release/example"
CLI_BIN="build/gmake_clang/bin/x86_64/Release/xatlas_cli"
MODELS_DIR="models"
VERBOSE=""

# Check if example binary exists
if [ ! -f "$EXAMPLE_BIN" ]; then
    echo "Error: Example binary not found at $EXAMPLE_BIN"
    echo "Please build the project first using: ./bin/premake.sh && cd build/gmake_clang && make"
    exit 1
fi

# Parse arguments
MODEL_NAME=""
USE_CLI=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --cli)
            USE_CLI=true
            shift
            ;;
        -verbose|--verbose)
            VERBOSE="--verbose"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [model_name] [options]"
            echo ""
            echo "Available models in $MODELS_DIR:"
            ls -1 "$MODELS_DIR"/*.obj 2>/dev/null | sed 's|.*/||' | sed 's|\.obj||' || echo "No .obj files found"
            echo ""
            echo "Options:"
            echo "  --cli       Use enhanced CLI tool with full options"
            echo "  -verbose    Enable verbose output"
            echo "  -h, --help  Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 cube"
            echo "  $0 gazebo -verbose"
            echo "  $0 Gargoyle_quadriflow"
            echo "  $0 cube --cli"
            echo ""
            echo "For full CLI options, run: $CLI_BIN --help"
            exit 0
            ;;
        -*)
            echo "Error: Unknown option '$1'"
            echo "Use '$0 --help' for usage information"
            exit 1
            ;;
        *)
            if [ -z "$MODEL_NAME" ]; then
                MODEL_NAME="$1"
            else
                echo "Error: Multiple model names specified: '$MODEL_NAME' and '$1'"
                exit 1
            fi
            shift
            ;;
    esac
done

# Check if model name was provided
if [ -z "$MODEL_NAME" ]; then
    echo "Error: No model name provided"
    echo "Usage: $0 [model_name] [options]"
    echo "Use '$0 --help' for more information"
    exit 1
fi

# Construct full path to model
MODEL_PATH="$MODELS_DIR/${MODEL_NAME}.obj"

# Check if model file exists
if [ ! -f "$MODEL_PATH" ]; then
    echo "Error: Model file not found: $MODEL_PATH"
    echo "Available models:"
    ls -1 "$MODELS_DIR"/*.obj 2>/dev/null | sed 's|.*/||' | sed 's|\.obj||' || echo "No .obj files found"
    exit 1
fi

# Run the example
if [ "$USE_CLI" = true ]; then
    echo "Running xatlas CLI on: $MODEL_PATH"
    echo "Command: $CLI_BIN $MODEL_PATH $VERBOSE"
    echo ""
    
    cd build/gmake_clang
    mkdir -p ../../tests/output
    ./bin/x86_64/Release/xatlas_cli "../../$MODEL_PATH" $VERBOSE
    # Move output files to tests/output directory
    mv output.obj ../../tests/output/ 2>/dev/null || true
    mv atlas_*.tga ../../tests/output/ 2>/dev/null || true
else
    echo "Running xatlas on: $MODEL_PATH"
    echo "Command: $EXAMPLE_BIN $MODEL_PATH $VERBOSE"
    echo ""
    
    cd build/gmake_clang
    mkdir -p ../../tests/output
    ./bin/x86_64/Release/example "../../$MODEL_PATH" $VERBOSE
    # Move output files to tests/output directory
    mv example_output.obj ../../tests/output/ 2>/dev/null || true
    mv example_*.tga ../../tests/output/ 2>/dev/null || true
fi

# Check if successful
if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Success! Output files generated:"
    if [ "$USE_CLI" = true ]; then
        ls -la output.obj atlas_*.tga 2>/dev/null || echo "No output files found"
    else
        ls -la example_output.obj example_*.tga 2>/dev/null || echo "No output files found"
    fi
else
    echo ""
    echo "❌ Error: xatlas processing failed"
    exit 1
fi
