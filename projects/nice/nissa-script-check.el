;;; nissa-script-check.el --- NICE syntax diagnostics -*- lexical-binding: t; -*-

(require 'flymake)
(require 'nissa-script-completion)

(defvar-local nissa-script-check-process nil
  "Syntax checker currently running for this buffer.")

(defvar-local nissa-script-check-request nil
  "Pending syntax request in the shared persistent worker.")

(defvar-local nissa-script-check-generation 0
  "Generation used to discard canceled or superseded syntax responses.")

(defun nissa-script-check-cancel ()
  "Cancel this buffer's previous syntax check."
  (cl-incf nissa-script-check-generation)
  (nissa-script-session-cancel nissa-script-check-request)
  (setq nissa-script-check-request nil)
  (when (processp nissa-script-check-process)
    (let ((process nissa-script-check-process))
      (setq nissa-script-check-process nil)
      (when (process-live-p process)
        (delete-process process)))))

(defun nissa-script-check-make-diagnostic (source begin end line column detail origin)
  "Map SOURCE's byte range BEGIN..END to a diagnostic for ORIGIN.
LINE, COLUMN, and DETAIL identify errors in preceding context files."
  (save-excursion
    (save-restriction
      (widen)
      (let ((position (point-min))
            (finish (point-min))
            (message detail))
        (if (equal source origin)
            (let ((bytes (encode-coding-string
                          (buffer-substring-no-properties (point-min) (point-max)) 'utf-8 t)))
              (setq position (+ (point-min) (length (decode-coding-string
                                                     (substring bytes 0 (min begin (length bytes))) 'utf-8-unix t))))
              (setq finish (+ (point-min) (length (decode-coding-string
                                                   (substring bytes 0 (min end (length bytes))) 'utf-8-unix t)))))
          (setq message (format "Context %s:%s:%s: %s" source line column detail)))
        ;; Anchor EOF errors to the last character so Flymake can underline them
        (when (and (= position (point-max)) (> position (point-min)))
          (setq position (1- position)))
        (flymake-make-diagnostic
         (current-buffer) position (min (point-max) (max finish (1+ position))) :error
         (if (string-empty-p message) "NICE syntax checker failed without a diagnostic" message))))))

(defun nissa-script-check-diagnostics (text origin)
  "Decode all ranged checker diagnostics from TEXT for ORIGIN.
The response is Lisp data, never evaluated.  Legacy text errors remain readable."
  (let ((data (condition-case nil (car (read-from-string text)) (error :invalid))))
    (if (and (listp data)
             (cl-every (lambda (entry)
                         (pcase entry
                           (`(,source ,begin ,end ,line ,column ,message)
                            (and (stringp source) (natnump begin) (natnump end) (<= begin end)
                                 (natnump line) (natnump column) (stringp message)))))
                       data))
        (mapcar (lambda (entry)
                  (apply #'nissa-script-check-make-diagnostic (append entry (list origin))))
                data)
      (if (string-match
           "Error: \\(.*\\):\\([0-9]+\\):\\([0-9]+\\) (offset \\([0-9]+\\)): \\(.*\\)" text)
          (let ((source (match-string 1 text))
                (line (string-to-number (match-string 2 text)))
                (column (string-to-number (match-string 3 text)))
                (offset (string-to-number (match-string 4 text)))
                (message (match-string 5 text)))
            (list (nissa-script-check-make-diagnostic source offset offset line column message origin)))
        (list (nissa-script-check-make-diagnostic origin 0 0 1 1 (string-trim text) origin))))))

(defun nissa-script-check-flymake (report-fn &rest _ignored)
  "Check the buffer through its shared worker, or the optional one-shot transport."
  (if nissa-script-use-persistent-session
      (nissa-script-check-flymake-session report-fn)
    (nissa-script-check-flymake-once report-fn)))

(defun nissa-script-check-flymake-session (report-fn)
  "Publish diagnostics from the retained AST without executing script code."
  (nissa-script-check-cancel)
  (let ((source (current-buffer))
        (generation nissa-script-check-generation)
        (tick (buffer-chars-modified-tick))
        (context (copy-tree nissa-script-completion-context))
        (directory default-directory)
        (origin (or buffer-file-name "<script>")))
    (condition-case error-data
        (setq nissa-script-check-request
              (nissa-script-session-request
               "check" (save-restriction (widen) (buffer-string)) 0
               (lambda (status payload message)
                 (when (buffer-live-p source)
                   (with-current-buffer source
                     (when (= generation nissa-script-check-generation)
                       (setq nissa-script-check-request nil))
                     (when (and (= generation nissa-script-check-generation)
                                (= tick (buffer-chars-modified-tick))
                                (equal context nissa-script-completion-context)
                                (equal directory default-directory)
                                (equal origin (or buffer-file-name "<script>")))
                       (setq nissa-script-check-request nil)
                       (funcall report-fn
                                (nissa-script-check-diagnostics
                                 (if (zerop status) (prin1-to-string payload) message) origin))))))))
      (error
       (funcall report-fn (nissa-script-check-diagnostics (error-message-string error-data) origin))))))

(defun nissa-script-check-flymake-once (report-fn)
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
               (command (append (list nissa-script-completion-program "--check-syntax" "--check-format" "emacs"
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
    (nissa-script-check-cancel)
    (unless (bound-and-true-p nissa-script-completion-mode)
      (nissa-script-session-stop))))

(provide 'nissa-script-check)
;;; nissa-script-check.el ends here
