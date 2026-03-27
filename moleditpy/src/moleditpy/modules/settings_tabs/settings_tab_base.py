from PyQt6.QtWidgets import QWidget

class SettingsTabBase(QWidget):
    """Base class for all settings tabs."""
    def __init__(self, default_settings, parent=None):
        super().__init__(parent)
        self.default_settings = default_settings

    def update_ui(self, settings_dict):
        """Update UI based on settings dictionary. Must be implemented by subclass."""
        raise NotImplementedError

    def get_settings(self):
        """Get settings values from the current UI. Must be implemented by subclass."""
        raise NotImplementedError

    def reset_to_defaults(self):
        """Reset only the settings of the current tab to defaults."""
        self.update_ui(self.default_settings)
