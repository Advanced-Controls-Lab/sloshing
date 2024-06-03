% net = importNetworkFromONNX("onnx_model_14.onnx");
% [240,320,3]
% cam=webcam(1);
vid = VideoWriter('vid','MPEG-4');
vid.FrameRate = 30;
open(vid);

tic
for i=1:500
    og_img=snapshot(cam);

    og_img=imresize(og_img,[240,320]);
    og_img=double(og_img)/255.0;
    
    output=predict(net,og_img);
    
    og_img=permute(og_img,[3,1,2]);
    a=og_img(1,:,:)+output;
    og_img(1,:,:)=a./max(a,[],'all');
    
    og_img=permute(og_img,[2,3,1]);
    im = uint8(round(og_img*255));
    % imshow(im)
    writeVideo(vid,im2frame(im));
end
toc;
close(vid);
