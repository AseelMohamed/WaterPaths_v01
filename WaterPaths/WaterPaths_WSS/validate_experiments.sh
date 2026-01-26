#!/bin/bash
# Quick validation script to test all experiment modes
# Run this after compilation to verify the experiment framework works correctly

echo "========================================="
echo "WaterPaths Experiment Framework Validator"
echo "========================================="
echo ""

# Check if executable exists
if [ ! -f "./WaterPaths_WSS" ]; then
    echo "ERROR: WaterPaths_WSS executable not found. Please compile first."
    exit 1
fi

# Check command line help includes -E flag
echo "1. Checking if -E flag is documented in help..."
./WaterPaths_WSS --help 2>&1 | grep -q "\-E: Experiment mode"
if [ $? -eq 0 ]; then
    echo "   ✓ -E flag found in help text"
else
    echo "   ✗ WARNING: -E flag not found in help text"
fi
echo ""

# Test invalid experiment mode
echo "2. Testing invalid experiment mode..."
./WaterPaths_WSS -E 5 -d InputFiles/ 2>&1 | grep -q "Invalid experiment mode"
if [ $? -eq 0 ]; then
    echo "   ✓ Invalid mode correctly rejected"
else
    echo "   ✗ WARNING: Invalid mode not properly handled"
fi
echo ""

# Test each valid experiment mode (dry run without actual simulation)
echo "3. Testing valid experiment modes..."
for exp in 1 2 3 4; do
    echo "   Testing Experiment $exp..."
    # This will fail if files don't exist, but that's okay - we just want to check parsing
    ./WaterPaths_WSS -E $exp -d InputFiles/ 2>&1 | head -1
done
echo ""

echo "========================================="
echo "Basic validation complete!"
echo ""
echo "Next steps:"
echo "1. Prepare a test solution file and input data"
echo "2. Run: ./WaterPaths_WSS -E 1 -m 0 -s test_solution.csv -d InputFiles/"
echo "3. Repeat for experiments 2, 3, and 4"
echo "4. Compare objective counts:"
echo "   - Experiments 1 & 2 should have 5 objectives"
echo "   - Experiments 3 & 4 should have 6 objectives"
echo "========================================="
