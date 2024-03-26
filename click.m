function click(~,eventData,handle_scatter,replace_mode)
% click is used as a callback ButtonDownFunction when selecting points on a surface
% click detect the click location and plot a pointer (defined by the
% handle_scatter) the selected point. replace_mode defines if the new
% pointer should be added to the old one or replace them.
Pt = eventData.IntersectionPoint;
hold on;
if replace_mode
    handle_scatter.XData=Pt(1);
    handle_scatter.YData=Pt(2);
    handle_scatter.ZData=Pt(3);
else
    handle_scatter.XData=[handle_scatter.XData Pt(1)];
    handle_scatter.YData=[handle_scatter.YData Pt(2)];
    handle_scatter.ZData=[handle_scatter.ZData Pt(3)];
end
end