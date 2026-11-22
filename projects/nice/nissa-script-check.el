;;; nissa-script-check.el --- NICE syntax diagnostics -*- lexical-binding: t; -*-

(require 'flymake)
(require 'nissa-script-completion)

(defvar-local nissa-script-check-process nil
  "Syntax checker currently running for this buffer.")

(defun nissa-script-check-cancel ()
  "Cancel this buffer's previous syntax check."
  (when (processp nissa-script-check-process)
    (let ((process nissa-script-check-process))
      (setq nissa-script-check-process nil)
      (when (process-live-p process)
        (delete-process process)))))

(defun nissa-script-check-diagnostics (text origin)
  "Convert checker TEXT to Flymake diagnostics for the current buffer and ORIGIN."
  (save-excursion
    (save-restriction
      (widen)
      (let ((position (point-min))
            (message (string-trim text)))
        (when (string-match
               "Error: \\(.*\\):\\([0-9]+\\):\\([0-9]+\\) (offset \\([0-9]+\\)): \\(.*\\)" text)
          (let ((source (match-string 1 text))
                (line (match-string 2 text))
                (column (match-string 3 text))
                (offset (string-to-number (match-string 4 text)))
                (detail (match-string 5 text)))
            (if (equal source origin)
                (let* ((bytes (encode-coding-string
                               (buffer-substring-no-properties (point-min) (point-max)) 'utf-8 t))
                       (prefix (substring bytes 0 (min offset (length bytes)))))
                  (setq position (+ (point-min) (length (decode-coding-string prefix 'utf-8-unix t))))
                  (setq message detail))
              (setq message (format "Context %s:%s:%s: %s" source line column detail)))))
        ;; Anchor EOF errors to the last character so Flymake can underline them
        (when (and (= position (point-max)) (> position (point-min)))
          (setq position (1- position)))
        (list (flymake-make-diagnostic
               (current-buffer) position (min (point-max) (1+ position)) :error
               (if (string-empty-p message) "NICE syntax checker failed without a diagnostic" message)))))))

(defun nissa-script-check-flymake (report-fn &rest _ignored)
  "Check the unsaved buffer asynchronously and send diagnostics to REPORT-FN."
  (nissa-script-check-cancel)
  (let ((source (current-buffer))
        (tick (buffer-chars-modified-tick))
        (origin (or buffer-file-name "<script>"))
        (output (generate-new-buffer " *NICE syntax check*")))
    (condition-case error-data
        (let* ((text (save-restriction
                       (widen)
                       (buffer-substring-no-properties (point-min) (point-max))))
               (command (append (list nissa-script-completion-program "--check-syntax"
                                      "--check-source" origin)
                                (nissa-script-completion-context-arguments)))
               (process
                (make-process
                 :name "NICE syntax check" :buffer output :stderr nil
                 :command command :connection-type 'pipe :noquery t
                 :coding '(utf-8-unix . utf-8-unix)
                 :sentinel
                 (lambda (process _event)
                   (when (memq (process-status process) '(exit signal))
                     (unwind-protect
                         (when (buffer-live-p source)
                           (with-current-buffer source
                             (when (eq process nissa-script-check-process)
                               (setq nissa-script-check-process nil)
                               (when (= tick (buffer-chars-modified-tick))
                                 (funcall report-fn
                                          (unless (zerop (process-exit-status process))
                                            (nissa-script-check-diagnostics
                                             (with-current-buffer output (buffer-string)) origin)))))))
                       (when (buffer-live-p output) (kill-buffer output))))))))
          (setq nissa-script-check-process process)
          (process-send-string process text)
          (process-send-eof process))
      (error
       (nissa-script-check-cancel)
       (when (buffer-live-p output) (kill-buffer output))
       (funcall report-fn
                (nissa-script-check-diagnostics (error-message-string error-data) origin))))))

(define-minor-mode nissa-script-check-mode
  "Underline NICE grammar errors in the unsaved buffer through Flymake."
  :lighter ""
  (if nissa-script-check-mode
      (progn
        (add-hook 'flymake-diagnostic-functions #'nissa-script-check-flymake nil t)
        (add-hook 'kill-buffer-hook #'nissa-script-check-cancel nil t)
        (flymake-mode 1))
    (remove-hook 'flymake-diagnostic-functions #'nissa-script-check-flymake t)
    (remove-hook 'kill-buffer-hook #'nissa-script-check-cancel t)
    (nissa-script-check-cancel)))

(provide 'nissa-script-check)
;;; nissa-script-check.el ends here
