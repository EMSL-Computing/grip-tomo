#!/bin/bash
# Reconstruct aligned tilt series into volumes using IMOD utilities.

tilt_bin="${IMOD_TILT:-tilt}"
trimvol_bin="${IMOD_TRIMVOL:-trimvol}"

for file in *.mrcs; do
	echo "$file"

	# Extract square dimension from IMOD header output.
	dims=$(header "$file" | grep "Number of columns, rows, sections" | awk '{print $7}')
	if [[ -z "$dims" ]]; then
		echo "Unable to determine image size for $file; skipping."
		continue
	fi
	echo "image size is: $dims x $dims"

	name=$(basename "$file")
	echo "$name"
	echo "Backprojecting $name"

	# Reconstruct the volume via IMOD tilt.
	"$tilt_bin" \
		-input "$file" \
		-output "$name.rec" \
		-IMAGEBINNED 1 \
		-TILTFILE *.rawtlt \
		-THICKNESS "$dims" \
		-RADIAL 0.35,0.035 \
		-FalloffIsTrueSigma 1 \
		-XAXISTILT 0.0 \
		-PERPENDICULAR \
		-MODE 2 \
		-FULLIMAGE "$dims","$dims" \
		-ActionIfGPUFails 1,2

	# Rotate the reconstruction by -90 degrees around the x-axis.
	"$trimvol_bin" "$name.rec" "$name.rec" -rx
done

rm *.rec~
