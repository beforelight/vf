[mag, phase, w] = bode(tfx1);

agl = phase(:) / 180 * pi;
mg = mag(:);
ww = w(:);
HH = mg .* cos(agl) + 1j * mg .* sin(agl);


HH = reshape(HH, [1, 1, max(size(HH))]);
Options.MaxIterations = 15;
Model = FastVF(ww, HH, 5, Options);

Hmodel = ComputeModelResponse(ww, Model.R0, Model.Rr, Model.Rc, Model.pr, Model.pc);
q = 1;
m = 1;
figure(33);
loglog(ww, abs(squeeze(Hmodel(q, m, :))), 'r-.', 'LineWidth', 1.5);
hold on;
grid on;
loglog(w(:), mag(:), 'b');
xlabel('\omega [rad/s]');
ylabel('Magnitude');
legend('Model', 'Samples');

figure(34);
semilogx(ww, angle_extern(angle(squeeze(Hmodel(q, m, :)))*180/pi), 'r-.', 'LineWidth', 1.5);
hold on;
grid on;
semilogx(w(:), agl(:)*180/pi, 'b');
ylabel('Phase [deg]');
xlabel('\omega [rad/s]');
legend('Model', 'Samples');


function agl = angle_extern(agl)
% 注意单位是deg
N = max(size(agl));
last = agl(1);
for i = 2:N
    while (last - 180 > agl(i)) || (agl(i) > last + 180)
        if agl(i) > last + 180
            agl(i) = agl(i) - 360;
        elseif last - 180 > agl(i)
            agl(i) = agl(i) + 360;
        end
    end
    last = agl(i);
end
end