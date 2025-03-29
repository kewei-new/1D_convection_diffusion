function print_progress(t, T, varargin)
% 打印动态进度条
% 参数：
%   t - 当前进度
%   T - 总进度
% 可选参数：
%   'BarLength' - 进度条长度（默认20）
%   'FillChar' - 填充字符（默认'█'）
%   'EmptyChar' - 空白字符（默认'-'）
%   'ShowPercent' - 显示百分比（默认true）
%   'Color' - 启用颜色（默认true）
%   'Prefix' - 前导文本（默认''）
%   'Suffix' - 后缀文本（默认''）

% 解析可选参数
p = inputParser;
addParameter(p, 'BarLength', 20, @isnumeric);
addParameter(p, 'FillChar', '█', @ischar);
addParameter(p, 'EmptyChar', '-', @ischar);
addParameter(p, 'ShowPercent', true, @islogical);
addParameter(p, 'Color', true, @islogical);
addParameter(p, 'Prefix', '', @ischar);
addParameter(p, 'Suffix', '', @ischar);
parse(p, varargin{:});

% 计算进度比例
progress = max(0, min(t / T, 1));
filled_len = round(p.Results.BarLength * progress);

% 构建进度条
bar_str = [repmat(p.Results.FillChar, 1, filled_len), ...
          repmat(p.Results.EmptyChar, 1, p.Results.BarLength - filled_len)];
text_str = sprintf('%s[%s]', p.Results.Prefix, bar_str);

% 添加百分比
if p.Results.ShowPercent
    text_str = [text_str, sprintf(' %.1f%%', progress*100)];
end

% 添加颜色（仅支持支持ANSI的终端）
% if p.Results.Color
%     if progress < 0.33
%         col_code = '\033[31m'; % 红色
%     elseif progress < 0.66
%         col_code = '\033[33m'; % 黄色
%     else
%         col_code = '\033[32m'; % 绿色
%     end
%     text_str = [col_code, text_str, '\033[0m'];
% end

% 添加后缀
text_str = [text_str, p.Results.Suffix];

% 打印进度（覆盖模式）
if t == 0
    fprintf('%s', text_str);
else
    fprintf(repmat('\b', 1, numel(text_str)+1)); % 删除前一行
    fprintf('%s', text_str);
end

% 完成后换行
if t >= T
    fprintf('\n');
end
end