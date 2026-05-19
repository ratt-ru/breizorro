"""Comprehensive tests for breizorro functionality"""

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from breizorro.breizorro import reproject_mask_to_reference
from breizorro.utils import apply_radial_cutoff, match_mask_shape


class TestMatchMaskShape:
    """Test cases for match_mask_shape function"""

    def test_crop_larger_mask_to_smaller_target(self):
        """Test cropping a larger mask to a smaller target shape"""
        large_mask = np.ones((100, 100))
        target_shape = (50, 50)
        result = match_mask_shape(large_mask, target_shape)

        assert result.shape == target_shape, "Shape mismatch!"
        assert np.array_equal(result, np.ones(target_shape)), "Values incorrect!"

    def test_extend_smaller_mask_to_larger_target(self):
        """Test extending a smaller mask to a larger target with zero-padding"""
        small_mask = np.ones((30, 30))
        target_shape = (80, 80)
        result = match_mask_shape(small_mask, target_shape)

        assert result.shape == target_shape, "Shape mismatch!"
        assert np.sum(result) == 30 * 30, "Extended mask should have zeros padding!"
        # Verify the original data is in the top-left corner
        assert np.array_equal(result[:30, :30], small_mask), "Data placement incorrect!"

    def test_same_shape_preservation(self):
        """Test that same-sized masks are preserved"""
        same_mask = np.ones((50, 50))
        target_shape = (50, 50)
        result = match_mask_shape(same_mask, target_shape)

        assert result.shape == target_shape
        assert np.array_equal(result, same_mask), "Should be identical!"

    def test_match_with_different_values(self):
        """Test matching with non-uniform mask values"""
        mask = np.arange(100).reshape(10, 10)
        target_shape = (15, 12)
        result = match_mask_shape(mask, target_shape)

        assert result.shape == target_shape
        assert np.array_equal(result[:10, :10], mask), "Original data should be preserved!"


class TestApplyRadialCutoff:
    """Test cases for apply_radial_cutoff function"""

    def test_radial_cutoff_basic(self):
        """Test basic radial cutoff application"""
        test_mask = np.ones((100, 100))
        radius = 30
        result = apply_radial_cutoff(test_mask, radius)

        assert result.shape == test_mask.shape, "Shape should be preserved!"

        # Check that center is still 1 and corners are 0
        center_y, center_x = 50, 50
        assert result[center_y, center_x] == 1.0, "Center should be 1!"
        assert result[0, 0] == 0.0, "Corner should be 0!"

        # Count non-zero pixels (should be approximately pi*r^2)
        non_zero_count = np.count_nonzero(result)
        expected_count = np.pi * radius**2
        assert abs(non_zero_count - expected_count) < 100, "Radial cutoff area incorrect!"

    def test_no_radial_cutoff_none(self):
        """Test that None radial cutoff returns original mask"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, None)

        assert np.array_equal(result, test_mask), "Should return original!"

    def test_no_radial_cutoff_zero_radius(self):
        """Test that zero radius returns original mask"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, 0)

        assert np.array_equal(result, test_mask), "Zero radius should return original!"

    def test_radial_cutoff_negative_radius(self):
        """Test that negative radius returns original mask"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, -1)

        assert np.array_equal(result, test_mask), "Negative radius should return original!"

    def test_large_radius_cutoff(self):
        """Test large radius cutoff covering most of the image"""
        test_mask = np.ones((100, 100))
        result = apply_radial_cutoff(test_mask, 100)

        total_ones = np.count_nonzero(result)
        assert total_ones > 9000, "Large radius should cover most of image!"

    def test_radial_cutoff_circular_pattern(self):
        """Test that radial cutoff creates a circular pattern"""
        test_mask = np.ones((101, 101))  # Odd size for exact center
        radius = 25
        result = apply_radial_cutoff(test_mask, radius)

        center_y, center_x = 50, 50

        # Points within radius should be non-zero
        assert result[center_y, center_x] == 1.0
        assert result[center_y + 10, center_x] == 1.0
        assert result[center_y, center_x + 10] == 1.0

        # Points outside radius should be zero
        assert result[0, 0] == 0.0
        assert result[100, 100] == 0.0

    def test_radial_cutoff_with_nonbinary_mask(self):
        """Test radial cutoff on non-binary mask (should preserve multiplication)"""
        test_mask = np.full((100, 100), 5.0)
        radius = 30
        result = apply_radial_cutoff(test_mask, radius)

        center_y, center_x = 50, 50
        assert result[center_y, center_x] == 5.0, "Center should preserve value!"
        assert result[0, 0] == 0.0, "Corner should be zeroed!"


class TestMaskOperations:
    """Test cases for mask merging and subtraction"""

    def test_mask_addition(self):
        """Test basic mask addition"""
        mask1 = np.zeros((50, 50))
        mask1[10:20, 10:20] = 1

        mask2 = np.zeros((50, 50))
        mask2[15:25, 15:25] = 1

        result = mask1 + mask2
        # Overlapping region should have value 2
        assert result[16, 16] == 2, "Overlapping region should sum!"
        # Non-overlapping should be 1
        assert result[11, 11] == 1, "Non-overlapping should be 1!"
        assert result[21, 21] == 1, "Non-overlapping should be 1!"

    def test_mask_subtraction(self):
        """Test mask subtraction"""
        mask1 = np.ones((50, 50))
        mask2 = np.zeros((50, 50))
        mask2[10:40, 10:40] = 1

        # Subtract mask2 from mask1
        mask1[mask2 != 0] = 0

        # Center should be zeroed
        assert mask1[25, 25] == 0, "Subtracted region should be 0!"
        # Edges should still be 1
        assert mask1[5, 5] == 1, "Non-subtracted region should be 1!"


class TestIslandOperations:
    """Test cases for island detection and manipulation"""

    def test_island_counting(self):
        """Test counting disconnected islands"""
        from scipy.ndimage import label

        mask = np.zeros((50, 50))
        # Create 3 separate islands
        mask[5:10, 5:10] = 1
        mask[20:25, 20:25] = 1
        mask[35:40, 35:40] = 1

        labeled, num_features = label(mask)

        assert num_features == 3, f"Should detect 3 islands, found {num_features}!"

    def test_island_labeling(self):
        """Test island labeling with unique IDs"""
        from scipy.ndimage import label

        mask = np.zeros((50, 50))
        mask[5:10, 5:10] = 1
        mask[20:25, 20:25] = 1

        labeled, num_features = label(mask)

        # Check that islands have different labels
        label1 = labeled[7, 7]
        label2 = labeled[22, 22]

        assert label1 != label2, "Islands should have different labels!"
        assert label1 > 0 and label2 > 0, "Labels should be positive!"

    def test_island_removal(self):
        """Test removing specific islands"""
        from scipy.ndimage import label

        mask = np.zeros((50, 50))
        mask[5:10, 5:10] = 1
        mask[20:25, 20:25] = 1
        mask[35:40, 35:40] = 1

        labeled, num_features = label(mask)

        # Remove island 2
        island_to_remove = labeled[22, 22]
        mask[labeled == island_to_remove] = 0

        # Re-count islands
        labeled_new, num_features_new = label(mask)
        assert num_features_new == 2, "Should have 2 islands after removal!"


class TestMorphologicalOperations:
    """Test cases for dilation and erosion"""

    def test_dilation(self):
        """Test mask dilation"""
        from scipy.ndimage import binary_dilation

        mask = np.zeros((50, 50))
        mask[24:26, 24:26] = 1  # Small 2x2 square

        # Dilate by 1 pixel
        dilated = binary_dilation(mask, iterations=1)

        # Dilated mask should be larger
        assert np.sum(dilated) > np.sum(mask), "Dilation should increase mask size!"

    def test_erosion(self):
        """Test mask erosion"""
        from scipy.ndimage import binary_erosion

        mask = np.zeros((50, 50))
        mask[10:40, 10:40] = 1  # Large square

        # Erode by 1 iteration
        eroded = binary_erosion(mask, iterations=1)

        # Eroded mask should be smaller
        assert np.sum(eroded) < np.sum(mask), "Erosion should decrease mask size!"

    def test_fill_holes(self):
        """Test filling holes in mask"""
        from scipy.ndimage import binary_fill_holes

        # Create a mask with a hole
        mask = np.ones((50, 50))
        mask[20:30, 20:30] = 0  # Hole in the center

        filled = binary_fill_holes(mask)

        # Hole should be filled
        assert filled[25, 25] == 1, "Hole should be filled!"
        assert np.sum(filled) == 50 * 50, "All pixels should be 1!"


class TestMinimumSize:
    """Test cases for minimum island size filtering"""

    def test_remove_small_islands(self):
        """Test removing islands smaller than threshold"""
        from scipy.ndimage import label
        from scipy.ndimage import sum as ndsum

        mask = np.zeros((100, 100))
        # Large island
        mask[10:40, 10:40] = 1
        # Small island
        mask[50:52, 50:52] = 1  # Only 4 pixels

        labeled, num_features = label(mask)

        # Filter by size
        minimum_size = 10
        island_areas = np.array(ndsum(mask, labeled, np.arange(labeled.max() + 1)))
        min_mask = island_areas >= minimum_size
        filtered_mask = min_mask[labeled.ravel()].reshape(labeled.shape)

        # Small island should be removed
        assert filtered_mask[51, 51] == 0, "Small island should be removed!"
        # Large island should remain
        assert filtered_mask[25, 25] == 1, "Large island should remain!"


class TestBinaryConversion:
    """Test cases for binary mask conversion"""

    def test_make_binary(self):
        """Test converting labeled mask to binary"""
        mask = np.array([[0, 1, 2], [3, 4, 0], [0, 5, 6]])

        binary_mask = (mask != 0).astype(int)

        # All non-zero values should become 1
        assert binary_mask[0, 1] == 1
        assert binary_mask[0, 2] == 1
        assert binary_mask[1, 0] == 1
        # Zeros should remain 0
        assert binary_mask[0, 0] == 0
        assert binary_mask[1, 2] == 0

    def test_mask_inversion(self):
        """Test mask inversion"""
        mask = np.array([[0, 1, 1], [0, 0, 1], [1, 1, 0]])

        inverted = (mask == 0).astype(int)

        # Zeros should become 1, ones should become 0
        assert inverted[0, 0] == 1
        assert inverted[0, 1] == 0
        assert inverted[2, 2] == 1


class TestWCSOperations:
    """Test cases for WCS coordinate operations"""

    def test_wcs_dropaxis(self):
        """Test dropping axes from WCS"""
        # Create a simple 4D WCS
        header = fits.Header()
        header["NAXIS"] = 4
        header["NAXIS1"] = 100
        header["NAXIS2"] = 100
        header["NAXIS3"] = 1
        header["NAXIS4"] = 1
        header["CRPIX1"] = 50
        header["CRPIX2"] = 50
        header["CRPIX3"] = 1
        header["CRPIX4"] = 1
        header["CRVAL1"] = 0.0
        header["CRVAL2"] = 0.0
        header["CRVAL3"] = 1.0
        header["CRVAL4"] = 1.0
        header["CDELT1"] = -0.001
        header["CDELT2"] = 0.001
        header["CDELT3"] = 1.0
        header["CDELT4"] = 1.0
        header["CTYPE1"] = "RA---SIN"
        header["CTYPE2"] = "DEC--SIN"
        header["CTYPE3"] = "FREQ"
        header["CTYPE4"] = "STOKES"

        wcs = WCS(header)

        # Drop to 2D
        while len(wcs.array_shape) > 2:
            wcs = wcs.dropaxis(len(wcs.array_shape) - 1)

        assert len(wcs.array_shape) == 2, "Should have 2D WCS!"
        assert wcs.array_shape == (100, 100), "Shape should be preserved!"


class TestReprojectMaskToReference:
    """Test cases for the extracted reprojection helper"""

    def _build_header(self, size, crpix):
        header = fits.Header()
        header["NAXIS"] = 2
        header["NAXIS1"] = size[1]
        header["NAXIS2"] = size[0]
        header["CRPIX1"] = crpix[1]
        header["CRPIX2"] = crpix[0]
        header["CRVAL1"] = 0.0
        header["CRVAL2"] = 0.0
        header["CDELT1"] = 1.0
        header["CDELT2"] = 1.0
        header["CTYPE1"] = "RA---TAN"
        header["CTYPE2"] = "DEC--TAN"
        return header

    def test_returns_input_when_shape_and_wcs_match(self):
        mask = np.zeros((10, 10), dtype=float)
        mask[4, 4] = 1.0
        header = self._build_header((10, 10), (5, 5))
        wcs_ref = WCS(header)

        result = reproject_mask_to_reference(mask, header, mask, wcs_ref)

        assert result.shape == mask.shape
        assert np.array_equal(result, mask)

    def test_reprojects_to_reference_shape(self):
        source = np.zeros((5, 5), dtype=float)
        source[2, 2] = 1.0
        source_header = self._build_header((5, 5), (3, 3))

        reference = np.zeros((11, 11), dtype=float)
        reference_header = self._build_header((11, 11), (6, 6))
        reference_wcs = WCS(reference_header)

        result = reproject_mask_to_reference(source, source_header, reference, reference_wcs)

        assert result.shape == reference.shape
        assert np.count_nonzero(result) >= 1
        assert result.sum() == pytest.approx(1.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
