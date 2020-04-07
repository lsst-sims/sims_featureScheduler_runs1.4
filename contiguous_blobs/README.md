Looking at ways we can alter the blob selection behavior. See if making blobs more contiguous helps with the overall performance.

## smooth_reward

Smooth the reward function with a Gaussian filter. In theory, this should get rid of little pesky local maxima.

## radius_constrict

Force all the fields to be within a given radius of the maximum.

## weighting_distance

Give weight to how close a field is to the maximum when selecting which fields to include in a blob.

