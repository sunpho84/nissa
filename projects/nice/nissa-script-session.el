;;; nissa-script-session.el --- Persistent NICE AST queries -*- lexical-binding: t; -*-

(require 'cl-lib)
(require 'subr-x)

(defvar nissa-script-completion-program)
(defvar nissa-script-completion-context)

(defcustom nissa-script-use-persistent-session t
  "Reuse one NICE worker and its AST per script buffer."
  :type 'boolean :group 'nissa-script)

(defcustom nissa-script-session-timeout 10
  "Maximum seconds to wait for a synchronous completion response."
  :type 'number :group 'nissa-script)

(defvar-local nissa-script-session-process nil
  "Worker shared by completion, documentation, and Flymake in this buffer.")

(defun nissa-script-session-stop ()
  "Stop this buffer's worker and release its retained AST."
  (interactive)
  (when (processp nissa-script-session-process)
    (let ((process nissa-script-session-process))
      (setq nissa-script-session-process nil)
      (process-put process 'closing t)
      (clrhash (process-get process 'pending))
      (when (process-live-p process) (delete-process process))
      (when-let* ((errors (process-get process 'errors))
                  ((buffer-live-p errors)))
        (kill-buffer errors)))))

(defun nissa-script-session--filter (process text)
  "Decode complete response frames from PROCESS's TEXT without evaluating data."
  (let ((pending (concat (or (process-get process 'partial) "") text)) callbacks)
    (while (string-match "\n" pending)
      (let* ((line (substring pending 0 (match-beginning 0)))
             (next (substring pending (match-end 0)))
             (reply (condition-case nil (read-from-string line) (error nil))))
        (setq pending next)
        (pcase (car reply)
          (`(,id ,status ,revision ,parses ,payload ,message)
           (when (and (natnump id) (memq status '(0 1)) (natnump revision)
                      (natnump parses) (stringp message)
                      (string-empty-p (string-trim (substring line (cdr reply)))))
             (process-put process 'revision revision)
             (process-put process 'parses parses)
             (when-let* ((callback (gethash id (process-get process 'pending))))
               (remhash id (process-get process 'pending))
               (push (list callback status payload message) callbacks)))))))
    ;; A callback can start another query; publish framing state before dispatch
    (process-put process 'partial pending)
    (dolist (callback (nreverse callbacks))
      (apply (car callback) (cdr callback)))))

(defun nissa-script-session--sentinel (process _event)
  "Fail outstanding requests if PROCESS exits unexpectedly."
  (when (memq (process-status process) '(exit signal))
    (unless (process-get process 'closing)
      (let ((message (when-let* ((errors (process-get process 'errors))
                                ((buffer-live-p errors)))
                       (with-current-buffer errors (string-trim (buffer-string))))))
        (when (string-empty-p (or message ""))
          (setq message "NICE editor worker exited before replying"))
        (let ((callbacks (hash-table-values (process-get process 'pending))))
          (clrhash (process-get process 'pending))
          (dolist (callback callbacks) (funcall callback 1 nil message)))))
    (when-let* ((errors (process-get process 'errors))
                ((buffer-live-p errors)))
      (kill-buffer errors))))

(defun nissa-script-session--worker ()
  "Return this buffer's live worker, restarting after executable/directory changes."
  (let* ((program (or (executable-find nissa-script-completion-program)
                      nissa-script-completion-program))
         (attributes (file-attributes (file-truename program)))
         (key (list nissa-script-completion-program program default-directory
                    (file-attribute-modification-time attributes)
                    (file-attribute-size attributes)
                    (file-attribute-inode-number attributes))))
    (unless (and (process-live-p nissa-script-session-process)
                 (equal key (process-get nissa-script-session-process 'key)))
      (nissa-script-session-stop)
      (let ((errors (generate-new-buffer " *NICE editor errors*")))
        (condition-case error-data
            (let ((process (make-process
                            :name "NICE editor" :buffer nil :stderr errors
                            :command (list nissa-script-completion-program "--editor-server")
                            :connection-type 'pipe :noquery t :coding '(utf-8-unix . utf-8-unix)
                            :filter #'nissa-script-session--filter
                            :sentinel #'nissa-script-session--sentinel)))
              (process-put process 'key key)
              (process-put process 'errors errors)
              (process-put process 'sequence 0)
              (process-put process 'pending (make-hash-table :test #'eql))
              (setq nissa-script-session-process process))
          (error (kill-buffer errors) (signal (car error-data) (cdr error-data)))))
      (add-hook 'kill-buffer-hook #'nissa-script-session-stop nil t)
      (add-hook 'change-major-mode-hook #'nissa-script-session-stop nil t))
    nissa-script-session-process))

(defun nissa-script-session--pieces (text)
  "Read ordered context pieces and append the current buffer's TEXT.
Reading files each request detects external edits without relying on timestamps."
  (let ((index 0))
    (append
     (mapcar
      (lambda (piece)
        (cl-incf index)
        (pcase piece
          (`(script . ,value) (cons (format "<context[%d]>" index) value))
          (`(file . ,path)
           (let ((name (expand-file-name path)))
             (unless (file-regular-p name) (error "Cannot read context file %s" name))
             (cons path (with-temp-buffer
                          (let ((coding-system-for-read 'utf-8-unix))
                            (insert-file-contents name))
                          (buffer-string)))))
          (_ (error "Unknown NICE completion context kind: %S" (car piece)))))
      nissa-script-completion-context)
     (list (cons (or buffer-file-name "<script>") text)))))

(defun nissa-script-session-request (operation text offset callback)
  "Send OPERATION, TEXT, and UTF-8 OFFSET.
Call CALLBACK with status, data, and message.
Return (PROCESS . ID), which can be canceled without destroying the worker."
  (let* ((pieces (nissa-script-session--pieces text))
         (process (nissa-script-session--worker))
         (id (1+ (process-get process 'sequence)))
         (frame (format "NICE/1 %d %s %d %d\n" id operation offset (length pieces))))
    (dolist (piece pieces)
      (setq frame (concat frame
                          (format "%d %d\n" (string-bytes (encode-coding-string (car piece) 'utf-8 t))
                                  (string-bytes (encode-coding-string (cdr piece) 'utf-8 t)))
                          (car piece) (cdr piece))))
    (process-put process 'sequence id)
    (puthash id callback (process-get process 'pending))
    (condition-case error-data
        (process-send-string process frame)
      (error
       (remhash id (process-get process 'pending))
       (signal (car error-data) (cdr error-data))))
    (cons process id)))

(defun nissa-script-session-cancel (request)
  "Drop REQUEST's callback while allowing its parse to populate the shared cache."
  (when (and (consp request) (processp (car request)))
    (remhash (cdr request) (process-get (car request) 'pending))))

(defun nissa-script-session-query (text offset)
  "Synchronously complete TEXT at OFFSET through this buffer's retained AST."
  (let (done status payload message)
    (let* ((request (nissa-script-session-request
                     "complete" text offset
                     (lambda (code data error) (setq done t status code payload data message error))))
           (deadline (+ (float-time) nissa-script-session-timeout)))
      (while (and (not done) (process-live-p (car request)) (< (float-time) deadline))
        (accept-process-output (car request) 0.05))
      (unless done
        (nissa-script-session-cancel request)
        (let ((callbacks (hash-table-values (process-get (car request) 'pending))))
          (clrhash (process-get (car request) 'pending))
          (dolist (callback callbacks)
            (funcall callback 1 nil "NICE editor query timed out")))
        (nissa-script-session-stop)
        (setq message "NICE editor query timed out"))
      (if (and done (zerop status)) payload
        (message "NICE completion failed: %s" message)
        nil))))

(defun nissa-script-session-status ()
  "Report this buffer's worker PID, source revision, and actual parse count."
  (interactive)
  (if (process-live-p nissa-script-session-process)
      (message "NICE worker %s: revision %s, parses %s"
               (process-id nissa-script-session-process)
               (process-get nissa-script-session-process 'revision)
               (process-get nissa-script-session-process 'parses))
    (message "No NICE editor worker for this buffer")))

(provide 'nissa-script-session)
;;; nissa-script-session.el ends here
