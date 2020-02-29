function [Ao, p] = superconv2(maps, K, p)
% SUPERCONV3 Perform convolution of a set of 3D "images" (4D array) with a
% set of 2D kernels (3D array), giving a 4D output. 
% For convolutional neural networks, the third dimension can be summed
% (e.g. sum(Ao,3)) and then the activation function can be applied.
% For best performance, store the parameters structure, p, which contains
% all of the setup to quickly invoke the convolution.

% Setup kernel and parameters and return them to caller for speed (reuse)
if isempty(p)
    tmp = K(1,1,1,1);
    if (isa(gather(tmp), 'single'))
        PRECISION = 'single';
    else
        PRECISION = 'double';
    end
    
    p.Nk =      int32(size(K)); % dimension of kernel
    p.Nkeven =  int32(not(mod([p.Nk(1) p.Nk(2)], 2)));
    p.Nkodd =   int32(mod([p.Nk(1) p.Nk(2)], 2));
    p.Nkh =     ([p.Nk(1) p.Nk(2)] - p.Nkodd)/2;
    p.Nm =      int32([size(maps,1) size(maps,2)]); % dimension of image
    
    % pre-allocate output memory
    p.result = gpuArray.zeros(p.Nm(1)-p.Nk(1)+1, p.Nm(2)-p.Nk(2)+1, p.Nk(3), PRECISION);
    
    % load and setup cuda kernel
    p.cudakernel = parallel.gpu.CUDAKernel('superconv2.ptx','superconv2.cu','superconv2');
    p.cudakernel.GridSize = [p.Nm(1)-p.Nk(1)+1, p.Nm(2)-p.Nk(2)+1, 1]; % dimensionality of output images
    p.cudakernel.ThreadBlockSize = [p.Nk(3), 1, 1];
    
    p.Nk = gpuArray(p.Nk);
    p.Nkeven = gpuArray(p.Nkeven);
    p.Nkodd = gpuArray(p.Nkodd);
    p.Nkh = gpuArray(p.Nkh);
    p.Nm = gpuArray(p.Nm);
end

Ao = feval(p.cudakernel, p.result, maps, K, p.Nm(1), p.Nm(2), p.Nk(1),...
           p.Nk(2), p.Nkh(1), p.Nkh(2), p.Nkeven(1), p.Nkeven(2));

end

