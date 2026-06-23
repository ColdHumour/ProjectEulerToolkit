@echo off
setlocal enabledelayedexpansion

:: 请在这里修改为您自己的命令前缀（示例为 echo）
set "prefix=jupyter nbconvert --to html --template classic"

:: 提示用户输入匹配文本
set /p "pattern=请输入要匹配的文本: "
if "%pattern%"=="" (
    echo 错误：未输入任何文本。
    exit /b 1
)

:: 搜索当前目录下文件名包含 pattern 的文件（排除文件夹）
:: dir /b /a-d 仅列出文件名，不显示额外信息；*%pattern%* 实现包含匹配
:: 使用 for /f 循环捕获输出，因为 dir 默认按名字升序排序，第一个即为所需文件
set "found="
for /f "delims=" %%i in ('dir /b /a-d "*%pattern%*" 2^>nul') do (
    set "found=%%i"
    goto :found_one
)

:found_one
if "%found%"=="" (
    echo 未找到任何文件名包含 "%pattern%" 的文件。
    exit /b 1
)

echo 找到的第一个文件: "%found%"

:: 拼接命令并执行（自动为文件名添加双引号，以支持空格等特殊字符）
%prefix% "%found%"

endlocal