import sys
import re
import threading
import os
import tempfile
import customtkinter as ctk
from tkinter import filedialog, messagebox

# ==========================================
# 全局主题设置
# ==========================================
ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")

class TextRedirector:
    def __init__(self, text_widget): 
        self.text_widget = text_widget
    def write(self, text):
        if "Ignoring fixed y limits" in text: return
        self.text_widget.after(0, self._append_text, text)
    def _append_text(self, text):
        self.text_widget.configure(state="normal")
        self.text_widget.insert("end", text)
        self.text_widget.see("end")
        self.text_widget.configure(state="disabled")
    def flush(self): pass

# ==========================================
# 主程序窗口类
# ==========================================
class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        # 1. 初始隐藏窗口
        self.withdraw() 
        
        # 2. 基础设置
        self.title("高通量阻抗及阿伦尼乌斯曲线自动化分析程序 v1.5")
        self.geometry("1000x780")
        self.minsize(900, 650)
        
        # 3. 设置图标
        try:
            base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
            icon_path = os.path.join(base_path, "icon.ico")
            if os.path.exists(icon_path):
                self.iconbitmap(icon_path)
        except: pass

        # 4. 初始化变量
        self.selected_folder = ""
        self.ratio_widgets = {}
        self.core_engine = None  
        self.is_paused = False
        self.pause_event = threading.Event()
        self.pause_event.set()
        self.frames = {}
        
        # 5. 构建 UI 框架
        self.build_welcome_page()
        self.build_main_page()

        # 在 mainloop 启动后 50ms 执行初始化序列
        self.after(50, self.startup_sequence)

    def startup_sequence(self):
        """主线程执行：初始化页面、释放闪屏、加载引擎"""
        self.show_frame("WelcomePage")
        
        # 在主线程加载数学引擎（防止 Altair 在子线程导入时崩溃）
        try:
            import core_logic
            self.core_engine = core_logic
        except Exception as e:
            # 修复 NameError: 提前捕获错误消息
            err_msg = str(e)
            self.after(0, lambda msg=err_msg: print(f"[执行异常]: {msg}"))

        # 打包环境处理
        if "NUITKA_ONEFILE_PARENT" in os.environ:
            try:
                import nuitka_splash_screen
                nuitka_splash_screen.close()
            except ImportError:
                p = os.path.join(tempfile.gettempdir(), f"onefile_{os.environ['NUITKA_ONEFILE_PARENT']}_splash_feedback.tmp")
                if os.path.exists(p):
                    try: os.unlink(p)
                    except: pass
            
            # 缓冲一下再显示主界面，防止重叠黑边
            self.after(250, self.show_main_window)
        else:
            self.show_main_window()

    def show_main_window(self):
        self.deiconify()
        self.focus_force()

    def show_frame(self, page_name):
        frame = self.frames[page_name]
        frame.grid(row=0, column=0, sticky="nsew")
        frame.tkraise()

    # ==========================================
    # UI 构建模块
    # ==========================================
    def build_welcome_page(self):
        frame = ctk.CTkFrame(self)
        self.grid_rowconfigure(0, weight=1); self.grid_columnconfigure(0, weight=1)
        self.frames["WelcomePage"] = frame
        frame.grid_rowconfigure(0, weight=1); frame.grid_rowconfigure(5, weight=1); frame.grid_columnconfigure(0, weight=1)

        title = ctk.CTkLabel(frame, text="高通量阻抗及\n阿伦尼乌斯曲线自动化分析程序 v1.5", font=ctk.CTkFont(size=36, weight="bold"))
        title.grid(row=1, column=0, pady=(0, 25))

        badge = ctk.CTkFrame(frame, fg_color="transparent")
        badge.grid(row=2, column=0, pady=(0, 30))
        ctk.CTkLabel(badge, text=" 👤 独立开发：刘畅 ", font=ctk.CTkFont(size=18, weight="bold"), text_color="white", fg_color="#1f538d", corner_radius=10).pack(side="left", padx=10)
        ctk.CTkLabel(badge, text=" 理论指导：陈胜洲 ", font=ctk.CTkFont(size=18, weight="bold"), text_color="white", fg_color="#333333", corner_radius=10).pack(side="left", padx=10)

        # --------------------------------
        instruction_text = (
            " 使用声明: 仅用于学术验证和同行评议\n\n"
            " 适用范围: 当前仅支持聚合物电解质阻塞电池或类似体系\n\n"
            " 数据格式: 支持 .txt 或 .csv 格式\n\n"
            " 数据文件命名规范提示(仅当需要生成阿伦尼乌斯曲线时遵守): \n\n"
            " 确保文件名遵循 A-B-C 格式 (数字必须为自然数，数字后可接英文) \n\n"
            "  A: 代表电解质种类 (Type)\n"
            "  B: 代表平行电池序号 (Sample ID)\n"
            "  C: 代表测试温度摄氏度 (Temperature)\n\n"
            "示例: 1-2-40.txt 代表 “1号电解质 - 2号平行电池 - 40摄氏度” "
        )
        
        text_box = ctk.CTkTextbox(frame, width=700, height=260, font=ctk.CTkFont(size=15), fg_color=("gray90", "gray16"))
        text_box.insert("0.0", instruction_text); text_box.configure(state="disabled")
        text_box.grid(row=4, column=0, pady=20)

        ctk.CTkButton(frame, text="我已阅读，开始设置", width=250, height=45, command=lambda: self.show_frame("MainPage")).grid(row=5, column=0, pady=20, sticky="n")

    def build_main_page(self):
        frame = ctk.CTkFrame(self)
        self.frames["MainPage"] = frame
        frame.grid_rowconfigure(1, weight=1); frame.grid_columnconfigure(0, weight=1); frame.grid_columnconfigure(1, weight=1)
        top = ctk.CTkFrame(frame, fg_color="transparent")
        top.grid(row=0, column=0, columnspan=2, sticky="new", padx=20, pady=20)
        top.grid_columnconfigure(0, weight=1); top.grid_columnconfigure(1, weight=1)

        left = ctk.CTkFrame(top)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))
        ctk.CTkLabel(left, text="📁 1. 全局数据与主程序设置", font=ctk.CTkFont(size=16, weight="bold")).pack(pady=10, padx=20, anchor="w")
        ctk.CTkButton(left, text="选择数据文件夹", command=self.select_folder).pack(pady=5, padx=20, anchor="w")
        self.lbl_folder = ctk.CTkLabel(left, text="未选择路径...", text_color="gray"); self.lbl_folder.pack(padx=20, anchor="w")
        self.eis_pdf_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(left, text="生成拟合 PDF 图表", variable=self.eis_pdf_var).pack(pady=10, padx=20, anchor="w")
        self.eis_dpi_combo = ctk.CTkComboBox(left, values=["150", "300", "600"], width=200); self.eis_dpi_combo.set("300"); self.eis_dpi_combo.pack(padx=20, pady=5, anchor="w")

        right = ctk.CTkFrame(top)
        right.grid(row=0, column=1, sticky="nsew", padx=(10, 0))
        self.sub_enable_var = ctk.BooleanVar(value=False)
        ctk.CTkSwitch(right, text="🔥 2. 阿伦尼乌斯分析程序", font=ctk.CTkFont(size=16, weight="bold"), 
                      variable=self.sub_enable_var, command=self.toggle_sub_program).pack(pady=10, padx=20, anchor="w")
        self.sub_settings_frame = ctk.CTkFrame(right, fg_color="transparent"); self.sub_settings_frame.pack(fill="both", expand=True, padx=20)
        
        dim_f = ctk.CTkFrame(self.sub_settings_frame, fg_color="transparent"); dim_f.pack(fill="x")
        self.entry_thick = ctk.CTkEntry(dim_f, width=80); self.entry_thick.insert(0, "185"); self.entry_thick.pack(side="left", padx=2)
        ctk.CTkLabel(dim_f, text="μm").pack(side="left")
        self.entry_diam = ctk.CTkEntry(dim_f, width=80); self.entry_diam.insert(0, "15.8"); self.entry_diam.pack(side="left", padx=(10, 2))
        ctk.CTkLabel(dim_f, text="mm").pack(side="left")

        self.scroll_ratios = ctk.CTkScrollableFrame(self.sub_settings_frame, height=120); self.scroll_ratios.pack(fill="x", pady=10)
        self.sub_pdf_var = ctk.BooleanVar(value=False); ctk.CTkSwitch(self.sub_settings_frame, text="生成子程序 PDF 结果", variable=self.sub_pdf_var).pack(anchor="w")
        self.sub_dpi_combo = ctk.CTkComboBox(self.sub_settings_frame, values=["150", "300", "600"], width=200); self.sub_dpi_combo.set("300"); self.sub_dpi_combo.pack(pady=5, anchor="w")
        self.toggle_sub_program()

        bot = ctk.CTkFrame(frame)
        bot.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=20, pady=(0, 20))
        bot.grid_rowconfigure(2, weight=1); bot.grid_columnconfigure(0, weight=1)
        ctrl = ctk.CTkFrame(bot, fg_color="transparent"); ctrl.grid(row=0, column=0, sticky="ew", padx=10, pady=10)
        self.progress_bar = ctk.CTkProgressBar(ctrl); self.progress_bar.set(0); self.progress_bar.pack(side="left", fill="x", expand=True, padx=(0, 10))
        self.btn_run = ctk.CTkButton(ctrl, text="🚀 开始运行", width=120, command=self.start_run_thread); self.btn_run.pack(side="left", padx=5)
        self.btn_pause = ctk.CTkButton(ctrl, text="⏸️ 暂停", width=100, state="disabled", command=self.toggle_pause); self.btn_pause.pack(side="left", padx=5)
        self.lbl_status = ctk.CTkLabel(bot, text="系统就绪"); self.lbl_status.grid(row=1, column=0, sticky="w", padx=15)
        self.console = ctk.CTkTextbox(bot, font=("Consolas", 12), fg_color="white", text_color="black"); self.console.grid(row=2, column=0, sticky="nsew", padx=10, pady=(0, 10))
        self.console.configure(state="disabled")
        sys.stdout = TextRedirector(self.console); sys.stderr = sys.stdout

    def toggle_sub_program(self):
        if self.sub_enable_var.get(): self.sub_settings_frame.pack(fill="both", expand=True, padx=20)
        else: self.sub_settings_frame.pack_forget()

    def select_folder(self):
        f = filedialog.askdirectory()
        if f: self.selected_folder = f; self.lbl_folder.configure(text=f); self.scan_folder_for_ratios()

    def scan_folder_for_ratios(self):
        for w in self.scroll_ratios.winfo_children(): w.destroy()
        self.ratio_widgets.clear()
        if not os.path.exists(self.selected_folder): return
        # 修正：支持 .csv 和 .txt
        all_files =[f for f in os.listdir(self.selected_folder) if f.lower().endswith(('.txt', '.csv'))]
        ratio_temps = {}
        pattern = re.compile(r'^(\d+).*?-.*?-(\d+)\.(txt|csv)$', re.IGNORECASE)
        for f in all_files:
            m = pattern.match(f)
            if m:
                r, t = int(m.group(1)), int(m.group(2))
                if r not in ratio_temps: ratio_temps[r] = set()
                ratio_temps[r].add(t)
        if not all_files:
            ctk.CTkLabel(self.scroll_ratios, text="文件夹内无数据文件", text_color="#E57373").pack(pady=10)
            return
        if not ratio_temps:
            ctk.CTkLabel(self.scroll_ratios, text=f"已识别 {len(all_files)} 个数据文件\n(命名不规范，子程序已禁用)", text_color="gray").pack(pady=10)
            return
        ctk.CTkLabel(self.scroll_ratios, text=f"--- 已匹配子程序序列 ({len(all_files)}个文件) ---", font=ctk.CTkFont(size=12)).pack(pady=5)
        for r in sorted(ratio_temps.keys()):
            row = ctk.CTkFrame(self.scroll_ratios, fg_color="transparent"); row.pack(fill="x", pady=2)
            ctk.CTkLabel(row, text=f"序列 {r} ({len(ratio_temps[r])}温度点)", width=130, anchor="w").pack(side="left", padx=5)
            mode = ctk.CTkOptionMenu(row, values=["自动推断", "手动输入"], width=100); mode.pack(side="left", padx=5)
            ent = ctk.CTkEntry(row, width=60, placeholder_text="℃", state="disabled"); ent.pack(side="left")
            mode.configure(command=lambda choice, e=ent: e.configure(state="normal" if choice=="手动输入" else "disabled"))
            self.ratio_widgets[r] = {"mode": mode, "entry": ent}

    def toggle_pause(self):
        if self.is_paused: self.is_paused = False; self.btn_pause.configure(text="⏸️ 暂停", fg_color="#C07F00"); self.pause_event.set()
        else: self.is_paused = True; self.btn_pause.configure(text="▶️ 继续", fg_color="green"); self.pause_event.clear()

    def update_progress(self, c, t, m):
        self.pause_event.wait()
        self.after(0, lambda: (self.progress_bar.set(c/t if t>0 else 0), self.lbl_status.configure(text=m)))

    def extract_dpi(self, s):
        m = re.search(r'\d+', str(s)); return int(m.group()) if m else 300

    def start_run_thread(self):
        if not self.selected_folder: messagebox.showerror("错误", "请选择文件夹"); return
        run_sub = self.sub_enable_var.get(); tm_config = {}
        if run_sub:
            try:
                float(self.entry_thick.get().strip()); float(self.entry_diam.get().strip())
                for r, ws in self.ratio_widgets.items():
                    if ws["mode"].get() == "手动输入":
                        val = ws["entry"].get().strip()
                        if not val: raise ValueError(f"序列 {r} 相变点为空")
                        tm_config[r] = float(val)
                    else: tm_config[r] = "Auto"
            except ValueError as e: messagebox.showerror("参数错误", str(e)); return
        self.btn_run.configure(state="disabled"); self.btn_pause.configure(state="normal")
        self.progress_bar.set(0); self.is_paused = False; self.pause_event.set()
        t = threading.Thread(target=self.run_process, args=(
            self.eis_pdf_var.get(), self.extract_dpi(self.eis_dpi_combo.get()),
            run_sub, self.entry_thick.get().strip(), self.entry_diam.get().strip(),
            self.sub_pdf_var.get(), self.extract_dpi(self.sub_dpi_combo.get()), tm_config
        ))
        t.daemon = True; t.start()

    def run_process(self, ep, ed, rs, th, di, sp, sd, tm):
        try:
            self.core_engine.run_main_program_entry(self.selected_folder, ep, ed, self.update_progress)
            if rs: self.core_engine.run_sub_program_entry(self.selected_folder, th, di, tm, sp, sd, self.update_progress)
            self.update_progress(1, 1, "✅ 分析完成！")
            messagebox.showinfo("完成", "任务结束")
        except Exception as e: print(f"错误: {str(e)}")
        finally: self.after(0, lambda: (self.btn_run.configure(state="normal"), self.btn_pause.configure(state="disabled")))

if __name__ == "__main__":
    app = App()
    app.mainloop()