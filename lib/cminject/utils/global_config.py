#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of CMInject
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# If you use this program for scientific work, you should correctly reference it; see LICENSE file for details.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

"""
Code for storing and distributing global configuration values, i.e., values that multiple objects may or may not be
interested in, for example the numerical time-step or the number of spatial dimensions.

Objects can explicitly ask for stored values for a certain configuration key, or subscribe to updates of the values
of a certain configuration key. They may raise exceptions when they notice an incompatibility of the new configuration
value with their own implementation, to prevent ill-defined simulations being run.
"""

from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum
from typing import Any, Union, List, Iterable, Callable, Dict, Set

__all__ = ['ConfigKey', 'ConfigSubscriber', 'GlobalConfig']


class Singleton:
    """
    A generic implementation of a singleton class, i.e., a class of which only one instance can be created.
    """
    def __init__(self, klass):
        self.klass = klass
        self.instance = None

    def __call__(self, *args, **kwargs):
        if self.instance is None:
            self.instance = self.klass(*args, **kwargs)
        return self.instance


class ConfigKey(Enum):
    """
    An enumeration defining the configuration keys that :class:`GlobalConfig` subscribers can subscribe to.
    """
    NUMBER_OF_DIMENSIONS = 1
    TIME_STEP = 2


class ConfigSubscriber(ABC):
    """
    A class that can subscribe to :class:`GlobalConfig`'s changes.
    """
    @abstractmethod
    def config_change(self, key: ConfigKey, value: Any) -> None:
        """
        Will be called whenever the value of any subscribed key changes.
        Will be called once at the time of subscribing, *IF* the value for the subscribed key(s) is not None.

        :param key: The ConfigKey that the change occurred for.
        :param value: The new value of the configuration value stored for the key ``key``.
        :return: Nothing (unused).
        """
        pass


@Singleton
class GlobalConfig:
    """
    A pub-sub implementation of a global configuration storage that any 'interested' objects can
      * ask for the current value once
      * subscribe to
      * change (globally)

    Defined as a :class:`Singleton`.
    Used for propagating global configuration like the number of spatial dimensions or the time step
    (see also :class:`ConfigKey`).
    """
    def __init__(self):
        self._subscribers: Dict[ConfigKey, Set[ConfigSubscriber]] = defaultdict(set)
        self._config: Dict[ConfigKey, Any] = {k: None for k in ConfigKey}

    def set(self, key: ConfigKey, value: Any) -> None:
        """
        Sets the value corresponding to ``key`` to a new value.

        :param key: The ConfigKey to change the value for.
        :param value: The new value.
        :raises: KeyError if ``key`` is not a valid key.
        """
        self._config[key] = value
        for subscriber in self._subscribers.get(key, []):
            subscriber.config_change(key, self._config[key])

    def get(self, key: ConfigKey) -> Any:
        """
        Gets the current value corresponding to ``key``.

        :param key: The ConfigKey to get the value for.
        :return: The value corresponding to ``key`.
        :raises: KeyError if ``key`` is not a valid key.
        """
        return self._config[key]

    @staticmethod
    def _call_for_all_pairs(subscriber: Union[ConfigSubscriber, Iterable[ConfigSubscriber]],
                            key: Union[ConfigKey, List[ConfigKey]],
                            fn: Callable[[ConfigSubscriber, ConfigKey], Any]):
        if not isinstance(subscriber, Iterable):
            subscriber = [subscriber]
        if not isinstance(key, list):
            key = [key]

        for key_ in key:
            for subscriber_ in subscriber:
                fn(subscriber_, key_)

    def _subscribe_one(self, subscriber: ConfigSubscriber, key: ConfigKey):
        self._subscribers[key].add(subscriber)
        # notify once at the moment of subscription, if value for this key has already been set
        if self._config[key] is not None:
            subscriber.config_change(key, self._config[key])

    def subscribe(self, subscriber: Union[ConfigSubscriber, Iterable[ConfigSubscriber]],
                  key: Union[ConfigKey, List[ConfigKey]]) -> None:
        """
        Subscribes one or more ConfigSubscriber instances to changes corresponding to one or more ConfigKeys.

        :param subscriber: A subscriber or a list of subscribers. Must be instances of ConfigSubscriber.
        :param key: A ConfigKey or a list of ConfigKeys.
        """
        self._call_for_all_pairs(subscriber, key, self._subscribe_one)

    def _unsubscribe_one(self, subscriber: ConfigSubscriber, key: ConfigKey):
        self._subscribers[key].remove(subscriber)

    def unsubscribe(self, subscriber: Union[ConfigSubscriber, Iterable[ConfigSubscriber]],
                    key: Union[ConfigKey, List[ConfigKey]]) -> None:
        """
        Unsubscribes one or more ConfigSubscriber instances from changes corresponding to one or more ConfigKeys.

        :param subscriber: A subscriber or a list of subscribers. Must be instances of ConfigSubscriber.
        :param key: A ConfigKey or a list of ConfigKeys.
        """
        self._call_for_all_pairs(subscriber, key, self._unsubscribe_one)

### Local Variables:
### fill-column: 100
### truncate-lines: t
### End:
