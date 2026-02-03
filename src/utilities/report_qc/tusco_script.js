function bindTuscoModalHandlers() {
  // Use event delegation for open-modal buttons (works with paginated tables)
  document.addEventListener('click', function (event) {
    if (event.target.classList.contains('open-modal')) {
      var modalId = event.target.getAttribute('data-modal');
      var el = document.getElementById(modalId);
      if (el) el.style.display = 'block';
    }
  });

  // Use event delegation for close-modal buttons
  document.addEventListener('click', function (event) {
    if (event.target.classList.contains('close-modal')) {
      var modalId = event.target.getAttribute('data-modal');
      var el = document.getElementById(modalId);
      if (el) el.style.display = 'none';
    }
  });

  // Close modal when clicking outside the modal content
  window.addEventListener('click', function (event) {
    document.querySelectorAll('.modal').forEach(function (modal) {
      if (event.target === modal) {
        modal.style.display = 'none';
      }
    });
  });
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', bindTuscoModalHandlers);
} else {
  bindTuscoModalHandlers();
}
